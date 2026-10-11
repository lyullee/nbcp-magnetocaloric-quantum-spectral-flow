$ErrorActionPreference = 'Stop'

$projectRoot = Split-Path -Parent $PSScriptRoot
$octaveBin = Join-Path $PSScriptRoot 'octave-11.3.0\octave-11.3.0-w64\mingw64\bin'
$gitExe = (Get-Command git -ErrorAction Stop).Source
$gitRoot = Split-Path (Split-Path $gitExe -Parent) -Parent
$gitBin = Join-Path $gitRoot 'mingw64\bin'

$octaveCurl = Join-Path $octaveBin 'libcurl-4.dll'
$gitCurl = Join-Path $gitBin 'libcurl-4.dll'
if (-not (Test-Path -LiteralPath $octaveCurl) -or -not (Test-Path -LiteralPath $gitCurl)) {
    throw 'Octave or Git for Windows libcurl was not found.'
}
if ((Get-Item -LiteralPath $octaveCurl).VersionInfo.FileVersion -ne
    (Get-Item -LiteralPath $gitCurl).VersionInfo.FileVersion) {
    throw 'Octave and signed Git libcurl versions differ; refusing the shim.'
}

$files = @('libcurl-4.dll','libbrotlidec.dll','libbrotlicommon.dll',
           'libpsl-5.dll','zlib1.dll')
foreach ($name in $files) {
    $source = Join-Path $gitBin $name
    if (-not (Test-Path -LiteralPath $source)) { throw "Missing $source" }
    if ((Get-AuthenticodeSignature -LiteralPath $source).Status -ne 'Valid') {
        throw "$source does not have a valid Authenticode signature."
    }
}

$backup = Join-Path $octaveBin 'libcurl-4.dll.octave-original'
if (-not (Test-Path -LiteralPath $backup)) {
    Copy-Item -LiteralPath $octaveCurl -Destination $backup
}
foreach ($name in $files) {
    Copy-Item -LiteralPath (Join-Path $gitBin $name) -Destination (Join-Path $octaveBin $name) -Force
}

$octave = Join-Path $octaveBin 'octave-cli.exe'
& $octave --no-gui --quiet --eval "disp(version)"
if ($LASTEXITCODE -ne 0) { throw "Octave smoke test failed with $LASTEXITCODE" }

Write-Host "Signed compatibility shim installed inside $octaveBin"
Write-Host 'Windows security policy was not changed.'
