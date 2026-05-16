param(
  [int]$Port = 8000,
  [switch]$NoServe
)

$ErrorActionPreference = "Stop"
$scriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$docsPath = Join-Path $scriptDir "docs"
$appFile = Join-Path $scriptDir "app.R"
$buildScript = Join-Path $scriptDir "build_shinylive.R"

if (-not (Test-Path $appFile)) {
  Write-Error ("app.R not found at {0}" -f $appFile)
}

if (-not (Test-Path $buildScript)) {
  Write-Error ("build_shinylive.R not found at {0}" -f $buildScript)
}

$rscript = Get-Command "Rscript" -ErrorAction SilentlyContinue
if (-not $rscript) {
  $fallbacks = @(
    "C:\Program Files\R\R-4.3.1\bin\x64\Rscript.exe",
    "C:\Program Files\R\R-4.3.1\bin\Rscript.exe"
  )
  $resolved = $fallbacks | Where-Object { Test-Path $_ } | Select-Object -First 1
  if (-not $resolved) {
    Write-Error "Rscript not found."
  }
  $rscript = $resolved
} else {
  $rscript = $rscript.Source
}

Write-Host ("Building Shinylive app from {0} to {1}" -f $appFile, $docsPath)
& $rscript $buildScript

if ($LASTEXITCODE -ne 0) {
  Write-Error ("Shinylive build failed with exit code {0}" -f $LASTEXITCODE)
}

if ($NoServe) {
  Write-Host ("Build complete. Static site updated at {0}" -f $docsPath)
  exit 0
}

Write-Host ("Serving {0} on http://localhost:{1}/" -f $docsPath, $Port)
Start-Process -FilePath "python" -ArgumentList @("-m", "http.server", $Port) -WorkingDirectory $docsPath
