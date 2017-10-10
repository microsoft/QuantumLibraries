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
# Unlike before, we manually specify since the order matters.
$libDirectory = Join-Path $PSScriptRoot "Microsoft.Quantum.Canon"
# FIXME: This is redundant with the sources in the csproj...!
$qflatSources = @(
    # Find and include the standard library.
    (Find-QflatCompiler | Split-Path -Resolve | % { Join-Path $_ ..\..\..\..\Library\standard.qb } | Resolve-Path),
    # Provide stubs for primitive operations.
    "Stubs.qb",

    # Provide definitions of the identity and nop.
    "Identity.qb",

    # Endianness.qb contains newtype declarations that are needed more broadly,
    # so we include it first.
    "Endianness.qb",

    "DataStructures/Stack.qb",

    # Similarly with OracleTypes.qb, save for that it depends on OperationPow.qb.
    "OperationPow.qb",
    "OracleTypes.qb",

    "ApplyToEach.qb",
    "ApplyToRange.qb",
    "Arithmetic.qb",
    "Bind.qb",

    "IterativePhaseEstimation.qb",
    "QFT.qb",
    "QuantumPhaseEstimation.qb",
    "ShiftOp.qb",
    "With.qb",

    "Paulis.qb",

    # # QECC
    "Qecc/Types.qb",
    "Qecc/Utils.qb",
    "Qecc/BitFlipCode.qb"
) | ForEach-Object {
    if (-not [System.IO.Path]::IsPathRooted($_)) {
        Join-Path $libDirectory $_
    } else {
        $_
    }
}

$qflatSources | ConvertFrom-Qflat -Verbose:$VerbosePreference
if ($LASTEXITCODE -eq 0) {
    Invoke-MsBuild (Resolve-Path .\QbLibs.sln) -ShowBuildOutputInCurrentWindow
}
