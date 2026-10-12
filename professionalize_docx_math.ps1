param(
  [Parameter(Mandatory=$true)][string]$InputDocx
)

$ErrorActionPreference = 'Stop'
$inputPath = (Resolve-Path -LiteralPath $InputDocx).Path
$word = $null
$document = $null
try {
  $word = New-Object -ComObject Word.Application
  $word.Visible = $false
  $word.DisplayAlerts = 0
  $document = $word.Documents.Open($inputPath, $false, $false)
  $count = $document.OMaths.Count
  if ($count -gt 0) {
    $document.OMaths.BuildUp()
  }
  $document.Save()
  [pscustomobject]@{File=$inputPath; Equations=$count; Status='PASS'}
} finally {
  if ($document) { $document.Close(0) }
  if ($word) { $word.Quit() }
  [GC]::Collect()
  [GC]::WaitForPendingFinalizers()
}
