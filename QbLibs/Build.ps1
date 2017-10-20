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
    "Stubs.qb",

    "Math/Types.qb",
    "Math/Constants.qb",

    "IterateThroughCartesianProduct.qb",

    "Combinators/ApplyToRange.qb",
    "Enumeration/Iter.qb",

    # Diagnostics
    "Asserts/AssertQubit.qb",
    "Asserts/AssertOperationsEqualReferenced.qb",
    "Asserts/AssertOperationsEqualInPlace.qb",

    # Provide definitions of the identity and nop.
    "Identity.qb"

    # Endianness.qb contains newtype declarations that are needed more broadly,
    # so we include it first.
    "Endianness.qb",

    "DataStructures/Stack.qb"

    # Similarly with OracleTypes.qb, save for that it depends on OperationPow.qb.
    "Combinators/OperationPow.qb",
    "PhaseEstimation/Types.qb",
    "Arithmetic.qb",
    "Bind.qb",

    "QFT.qb",
    "PhaseEstimation/Quantum.qb",
    "PhaseEstimation/Iterative.qb",
    # "AmplitudeAmplification.qb"
    "ShiftOp.qb",
    "Combinators/With.qb",

    "Paulis.qb",

    # QECC
    "Qecc/Types.qb",
    "Qecc/Utils.qb",
    "Qecc/BitFlipCode.qb",
    "Qecc/5QubitCode.qb",
    "Qecc/7QubitCode.qb"
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
