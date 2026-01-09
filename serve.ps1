param(
  [int]$Port = 8000
)

$ErrorActionPreference = "Stop"
$scriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$docsPath = Join-Path $scriptDir "docs"

if (-not (Test-Path $docsPath)) {
  Write-Error "docs folder not found. Run shinylive::export(appdir='.', destdir='docs') first."
}

Write-Host ("Serving {0} on http://localhost:{1}/" -f $docsPath, $Port)
Start-Process -FilePath "python" -ArgumentList @("-m", "http.server", $Port) -WorkingDirectory $docsPath
