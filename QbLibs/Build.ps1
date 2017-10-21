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
    # This switch is off for now, since the solution doesn't yet build.
    [switch]
    $BuildSolution = $false
)

Import-Module -Force .\Invoke-Qbc.psm1
if (-not (Get-Module -ListAvailable Invoke-MsBuild)) {
    Install-Module Invoke-MsBuild
}
Import-Module Invoke-MsBuild

# FIXME: we'll want to do the following to use Visual Studio support.

# # Get a list of what we need to compile.
# # We do so by examining the *.csproj inside, looking for <None>
# # elements referencing *.qb files.
# $libDirectory = Join-Path $PSScriptRoot "Microsoft.Quantum.Canon"
# # Import the csproj as an XML document.
# $csproj = [xml](
#     Join-Path $libDirectory "Microsoft.Quantum.Canon.csproj" `
#     | ForEach-Object { Get-Content $_ })
# $qflatSources = $csproj.Project.ItemGroup `
#     | ForEach-Object { $_.None } `
#     | Select-Object -ExpandProperty Include `
#     | Where-Object { $_.EndsWith(".qb") } `
#     | Where-Object { $_.StartsWith("Phase") } `
#     | ForEach-Object { Join-Path $libDirectory $_ }

# Get a list of what we need to compile.
# Unlike before, we manually specify since the order matters.
$libDirectory = Join-Path $PSScriptRoot "Microsoft.Quantum.Canon"
# Find and include the standard library.
$qflatSources = @(Find-QflatPrelude)

$qflatSources += @(
    # Provide stubs for primitive operations.
    "Stubs.qb"

    "Math/NativeStubs.qb"

    "Combinators/RestrictToSubregister.qb"
    "Combinators/ApplyToEach.qb"

    "DataStructures/Pairs.qb"
    
    "Enumeration/Iter.qb"
    "Enumeration/Trotter.qb"

    "Simulation/Types.qb"
    "Simulation/SimulationTechniques.qb"
    "Simulation/EvolutionSetPauli.qb"
    #"Simulation/EvolutionSetFermionic.qb"

    "Simulation/ExampleH2.qb"
    "Simulation/ExampleIsing.qb"

    #"Simulation/Schedule.qb"

    #"Enumeration/Iter2.qb"

) | ForEach-Object {
    Join-Path $libDirectory $_
}

$qflatSources | ConvertFrom-Qflat

if ($BuildSolution) {
    if (-not (Get-Module -ListAvailable Invoke-MsBuild)) {
        Install-Module Invoke-MsBuild
    }
    Import-Module Invoke-MsBuild
    if ($LASTEXITCODE -eq 0) {
        Invoke-MsBuild (Resolve-Path .\QbLibs.sln) -ShowBuildOutputInCurrentWindow
    }
}
