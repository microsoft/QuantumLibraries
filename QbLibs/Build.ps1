##
# .SYNOPSIS
#     ./Build.ps1: This PowerShell script builds all QbLibs to generated C#
#     sources.
#
# .EXAMPLE
#     PS> $Env:QFLAT_PATH = "<path to Compiler.exe>"
#     PS> ./Build.ps1
##
param(
    [switch]
    $ParseOnly = $false
)

Import-Module -Force .\Invoke-Qbc.psm1

# Get a list of what we need to compile.
# Unlike before, we manually specify since the order matters.
$libDirectory = Join-Path $PSScriptRoot "Microsoft.Quantum.Canon"
$qflatSources = @(
    # Provide stubs for primitive operations.
    "Stubs.qb",

    # Endianness.qb contains newtype declarations that are needed more broadly,
    # so we include it first.
    "Endianness.qb",

    # Similarly with OracleTypes.qb, save for that it depends on OperationPow.qb.
    "OperationPow.qb",
    "OracleTypes.qb",

    "ApplyToEach.qb",
    "ApplyToRange.qb",
    "Arithmetic.qb",
    "IterativePhaseEstimation.qb",
    # # "QFT.qb", # QFT commented out in lieu of merging in martinro/ branch.
    "QuantumPhaseEstimation.qb",
    "ShiftOp.qb",
    "With.qb",

    "Paulis.qb",
    
    # # QECC
    "Qecc/Types.qb",
    "Qecc/Utils.qb",
    "Qecc/BitFlipCode.qb"
) | ForEach-Object {
    Join-Path $libDirectory $_
}

$qflatSources | ConvertFrom-Qflat -ParseOnly:$ParseOnly -Batch -Target Canon.g.cs
