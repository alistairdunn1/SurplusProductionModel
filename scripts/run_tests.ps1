# PowerShell wrapper to run R tests for SurplusProductionModel
param()

# Try common locations for Rscript
$possible = @(
  "$env:ProgramFiles\R",
  "$env:ProgramFiles\R\R-4.4.1\bin",
  "$env:ProgramFiles\R\R-4.4.0\bin",
  "$env:ProgramFiles\R\R-4.3.3\bin",
  "$env:ProgramFiles\R\R-4.3.2\bin",
  "$env:LOCALAPPDATA\Programs\R",
  "$env:LOCALAPPDATA\Programs\R\R-4.4.1\bin",
  "$env:LOCALAPPDATA\Programs\R\R-4.4.0\bin"
) | Where-Object { $_ -and (Test-Path $_) }

$RscriptPath = $null
foreach ($p in $possible) {
  $candidate = Join-Path $p "Rscript.exe"
  if (Test-Path $candidate) { $RscriptPath = $candidate; break }
  $cand = Get-ChildItem -Path $p -Recurse -Filter Rscript.exe -ErrorAction SilentlyContinue | Select-Object -First 1 -ExpandProperty FullName
  if ($cand) { $RscriptPath = $cand; break }
}

if (-not $RscriptPath) {
  Write-Error "Could not find Rscript.exe. Ensure R is installed and Rscript is on PATH."
  exit 1
}

$script = Join-Path $PSScriptRoot "run_tests.R"
& "$RscriptPath" "$script"
exit $LASTEXITCODE
