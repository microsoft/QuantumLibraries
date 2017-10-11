##
# .SYNOPSIS
#     ./Build.ps1: This PowerShell script builds all QbLibs to generated C#
#     sources.
#
# .EXAMPLE
#     PS> $Env:QFLAT_PATH = "<path to Compiler.exe>"
#     PS> ./Build.ps1
##
[CmdletBinding()]
param(
)

Import-Module -Force .\Invoke-Qbc.psm1
if (-not (Get-Module -ListAvailable Invoke-MsBuild)) {
    Install-Module Invoke-MsBuild
}
Import-Module Invoke-MsBuild

# Get a list of what we need to compile.
# We do so by examining the *.csproj inside, looking for <None>
# elements referencing *.qb files.
$libDirectory = Join-Path $PSScriptRoot "Microsoft.Quantum.Canon"
# Import the csproj as an XML document.
$csproj = [xml](
    Join-Path $libDirectory "Microsoft.Quantum.Canon.csproj" `
    | ForEach-Object { Get-Content $_ })
$qflatSources = $csproj.Project.ItemGroup `
    | ForEach-Object { $_.None } `
    | Select-Object -ExpandProperty Include `
    | Where-Object { $_.EndsWith(".qb") } `
    | ForEach-Object { Join-Path $libDirectory $_ }

# Find and include the standard library.
# $qflatSources += @(
#      (Find-QflatCompiler `
#         | Split-Path -Resolve `
#         | ForEach-Object { Join-Path $_ ..\..\..\..\Library\standard.qb } `
#         | Resolve-Path
#     )
# )

# "Compiling Q♭ sources:`n$(
#     "`n" -join ($qflatSources | ForEach-Object { "`t$_" })
# )" | Write-Verbose

# "Compiling Q♭ sources:`n$(
#     $qflatSources
# )" | Write-Verbose

$qflatSources | ConvertFrom-Qflat -Verbose:$VerbosePreference
if ($LASTEXITCODE -eq 0) {
    Invoke-MsBuild (Resolve-Path .\QbLibs.sln) -ShowBuildOutputInCurrentWindow
}
