$ErrorActionPreference = 'Stop'
$root = Split-Path -Parent $MyInvocation.MyCommand.Path
Set-Location -LiteralPath $root
$python = (Get-Command python).Source
$results = Join-Path $root 'results'
$statusPath = Join-Path $results 'n18_multimodel_pipeline_status.json'
$denseJobsPath = Join-Path $results 'n18_multimodel_dense_jobs.json'

function Write-Status([string]$stage, [string]$status, [object]$details) {
    [ordered]@{
        stage = $stage
        status = $status
        updated_at = (Get-Date).ToString('o')
        details = $details
    } | ConvertTo-Json -Depth 8 | Set-Content -LiteralPath $statusPath -Encoding utf8
}

Write-Status 'dense-refine' 'waiting' @{
    note = 'Waiting for the five model-specific common-grid workers.'
}
$denseJobs = Get-Content -LiteralPath $denseJobsPath -Raw | ConvertFrom-Json
while ($true) {
    $running = @($denseJobs | Where-Object {
        Get-Process -Id $_.pid -ErrorAction SilentlyContinue
    })
    if ($running.Count -eq 0) { break }
    Start-Sleep -Seconds 20
}

$denseErrors = @()
foreach ($job in $denseJobs) {
    $safe = $job.model.Replace('-', '')
    $checkpointPath = Join-Path $results "n18_multimodel_refine_checkpoint_${safe}.json"
    $traceback = (Test-Path -LiteralPath $job.stderr) -and
        ((Get-Content -LiteralPath $job.stderr -Raw) -match 'Traceback')
    if (-not (Test-Path -LiteralPath $checkpointPath)) {
        $denseErrors += [pscustomobject]@{model=$job.model; reason='missing checkpoint'}
        continue
    }
    $checkpoint = Get-Content -LiteralPath $checkpointPath -Raw | ConvertFrom-Json
    if ($checkpoint.rows.Count -ne 51) {
        $denseErrors += [pscustomobject]@{model=$job.model; reason="checkpoint has $($checkpoint.rows.Count) rows"}
    } elseif ($traceback) {
        $denseErrors += [pscustomobject]@{model=$job.model; reason='traceback in stderr'}
    }
}
if ($denseErrors.Count -gt 0) {
    Write-Status 'dense-refine' 'failed' @{
        models = @($denseErrors.model)
        failures = $denseErrors
        note = 'At least one dense worker lacks a complete 51-row checkpoint or raised a traceback.'
    }
    exit 1
}

$models = @('M3', 'M4-00', 'M4-02', 'M4-03', 'M4-01')
$convergenceJobs = @()
foreach ($model in $models) {
    $safe = $model.Replace('-', '')
    $stdout = "results\n18_${safe}_convergence_stdout.log"
    $stderr = "results\n18_${safe}_convergence_stderr.log"
    $process = Start-Process -FilePath $python `
        -ArgumentList @('run_n18_multimodel_cover.py', '--stage', 'convergence', '--models', $model) `
        -WorkingDirectory $root -WindowStyle Hidden `
        -RedirectStandardOutput $stdout -RedirectStandardError $stderr -PassThru
    $process.PriorityClass = 'BelowNormal'
    $convergenceJobs += [pscustomobject]@{
        model = $model
        pid = $process.Id
        stdout = $stdout
        stderr = $stderr
    }
}
$convergenceJobs | ConvertTo-Json | Set-Content `
    -LiteralPath (Join-Path $results 'n18_multimodel_convergence_jobs.json') -Encoding utf8
Write-Status 'convergence' 'running' @{
    models = $models
    jobs = $convergenceJobs
}

while ($true) {
    $running = @($convergenceJobs | Where-Object {
        Get-Process -Id $_.pid -ErrorAction SilentlyContinue
    })
    if ($running.Count -eq 0) { break }
    Start-Sleep -Seconds 20
}

$convergenceErrors = @()
foreach ($job in $convergenceJobs) {
    $safe = $job.model.Replace('-', '')
    $checkpointPath = Join-Path $results "n18_multimodel_convergence_checkpoint_${safe}.json"
    $traceback = (Test-Path -LiteralPath $job.stderr) -and
        ((Get-Content -LiteralPath $job.stderr -Raw) -match 'Traceback')
    if (-not (Test-Path -LiteralPath $checkpointPath)) {
        $convergenceErrors += [pscustomobject]@{model=$job.model; reason='missing checkpoint'}
        continue
    }
    $checkpoint = Get-Content -LiteralPath $checkpointPath -Raw | ConvertFrom-Json
    if (-not $checkpoint.models.PSObject.Properties[$job.model]) {
        $convergenceErrors += [pscustomobject]@{model=$job.model; reason='model absent from checkpoint'}
    } elseif ($traceback) {
        $convergenceErrors += [pscustomobject]@{model=$job.model; reason='traceback in stderr'}
    }
}
if ($convergenceErrors.Count -gt 0) {
    Write-Status 'convergence' 'failed' @{
        models = @($convergenceErrors.model)
        failures = $convergenceErrors
        note = 'At least one convergence worker lacks its model checkpoint or raised a traceback.'
    }
    exit 1
}

Write-Status 'compile' 'running' @{}
& $python 'verify_n18_multimodel_plan.py' `
    1> 'results\n18_multimodel_plan_verify_stdout.log' `
    2> 'results\n18_multimodel_plan_verify_stderr.log'
if ($LASTEXITCODE -ne 0) {
    Write-Status 'plan-verify' 'failed' @{ exit_code = $LASTEXITCODE }
    exit $LASTEXITCODE
}

& $python 'compile_n18_multimodel_cover.py' `
    1> 'results\n18_multimodel_compile_stdout.log' `
    2> 'results\n18_multimodel_compile_stderr.log'
if ($LASTEXITCODE -ne 0) {
    Write-Status 'compile' 'failed' @{ exit_code = $LASTEXITCODE }
    exit $LASTEXITCODE
}

& $python 'verify_n18_multimodel_cover.py' `
    1> 'results\n18_multimodel_verify_stdout.log' `
    2> 'results\n18_multimodel_verify_stderr.log'
if ($LASTEXITCODE -ne 0) {
    Write-Status 'verify' 'failed' @{ exit_code = $LASTEXITCODE }
    exit $LASTEXITCODE
}

& $python 'stress_test_n18_multimodel_design.py' `
    1> 'results\n18_multimodel_design_stdout.log' `
    2> 'results\n18_multimodel_design_stderr.log'
if ($LASTEXITCODE -ne 0) {
    Write-Status 'design-propagation' 'failed' @{ exit_code = $LASTEXITCODE }
    exit $LASTEXITCODE
}

& $python 'verify_n18_multimodel_design.py' `
    1> 'results\n18_multimodel_design_verify_stdout.log' `
    2> 'results\n18_multimodel_design_verify_stderr.log'
if ($LASTEXITCODE -ne 0) {
    Write-Status 'design-verify' 'failed' @{ exit_code = $LASTEXITCODE }
    exit $LASTEXITCODE
}

Write-Status 'manuscript-integration' 'running' @{}
& $python 'integrate_n18_multimodel_manuscript.py' `
    1> 'results\n18_multimodel_manuscript_stdout.log' `
    2> 'results\n18_multimodel_manuscript_stderr.log'
if ($LASTEXITCODE -ne 0) {
    Write-Status 'manuscript-integration' 'failed' @{ exit_code = $LASTEXITCODE }
    exit $LASTEXITCODE
}

& $python 'verify_manuscript_claims.py' `
    1> 'results\n18_multimodel_claims_verify_stdout.log' `
    2> 'results\n18_multimodel_claims_verify_stderr.log'
if ($LASTEXITCODE -ne 0) {
    Write-Status 'manuscript-claims-verify' 'failed' @{ exit_code = $LASTEXITCODE }
    exit $LASTEXITCODE
}

& $python 'verify_cryogenics_submission.py' `
    1> 'results\n18_multimodel_submission_verify_stdout.log' `
    2> 'results\n18_multimodel_submission_verify_stderr.log'
if ($LASTEXITCODE -ne 0) {
    Write-Status 'submission-verify' 'failed' @{ exit_code = $LASTEXITCODE }
    exit $LASTEXITCODE
}

$testStdout = Join-Path $results 'n18_multimodel_targeted_tests_stdout.log'
$testStderr = Join-Path $results 'n18_multimodel_targeted_tests_stderr.log'
Set-Content -LiteralPath $testStdout -Value '' -Encoding utf8
Set-Content -LiteralPath $testStderr -Value '' -Encoding utf8
foreach ($testScript in @(
    'test_general_truncated_magnon_operator.py',
    'test_stacked_lattice.py'
)) {
    $safeTest = [IO.Path]::GetFileNameWithoutExtension($testScript)
    $tempStdout = Join-Path $results "${safeTest}_stdout.tmp"
    $tempStderr = Join-Path $results "${safeTest}_stderr.tmp"
    $testProcess = Start-Process -FilePath $python `
        -ArgumentList @($testScript) -WorkingDirectory $root -WindowStyle Hidden `
        -RedirectStandardOutput $tempStdout -RedirectStandardError $tempStderr `
        -Wait -PassThru
    Add-Content -LiteralPath $testStdout -Value "[$testScript]"
    if (Test-Path -LiteralPath $tempStdout) {
        Add-Content -LiteralPath $testStdout -Value (Get-Content -LiteralPath $tempStdout -Raw)
        Remove-Item -LiteralPath $tempStdout
    }
    if (Test-Path -LiteralPath $tempStderr) {
        Add-Content -LiteralPath $testStderr -Value (Get-Content -LiteralPath $tempStderr -Raw)
        Remove-Item -LiteralPath $tempStderr
    }
    if ($testProcess.ExitCode -ne 0) {
        Write-Status 'targeted-tests' 'failed' @{
            script = $testScript
            exit_code = $testProcess.ExitCode
        }
        exit $testProcess.ExitCode
    }
}

& $python 'make_bulk_manifest.py' `
    1> 'results\n18_multimodel_manifest_stdout.log' `
    2> 'results\n18_multimodel_manifest_stderr.log'
if ($LASTEXITCODE -ne 0) {
    Write-Status 'manifest' 'failed' @{ exit_code = $LASTEXITCODE }
    exit $LASTEXITCODE
}

Write-Status 'complete' 'passed' @{
    summary = 'results/n18_multimodel_cover_summary.json'
    verification = 'results/n18_multimodel_cover_verification.json'
    figure = 'figures/n18_multimodel_cover.png'
    design_summary = 'results/n18_multimodel_design_stress_summary.json'
    design_verification = 'results/n18_multimodel_design_stress_verification.json'
    design_figure = 'figures/n18_multimodel_design_stress.png'
    manuscript_integration = 'results/n18_multimodel_manuscript_integration.json'
    manuscript_verification = 'results/manuscript_claim_verification.json'
    submission_verification = 'results/cryogenics_submission_verification.json'
    manifest = 'results/bulk_path_manifest_20261011.json'
}
