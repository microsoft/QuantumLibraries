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

# Get a list of what we need to compile.
# Unlike before, we manually specify since the order matters.
$libDirectory = Join-Path $PSScriptRoot "Microsoft.Quantum.Canon"
# Find and include the standard library.
$qflatSources = @(Find-QflatPrelude)

$qflatSources += @(
    # Provide stubs for primitive operations.
    "Stubs.qb"

    "Math/Types.qb",
    "Math/Constants.qb",

    "IterateThroughCartesianProduct.qb",
    
    "Combinators/ApplyToEach.qb",
    "Combinators/ApplyToRange.qb",
    "Combinators/OperationPow.qb",
    "Combinators/With.qb",

    # Diagnostics
    "Asserts/AssertQubit.qb",
    "Asserts/AssertOperationsEqualReferenced.qb",
    "Asserts/AssertOperationsEqualInPlace.qb",
	"Asserts/ClassicalAsserts.qb"

    # # Provide definitions of the identity and nop.
    "Identity.qb"

    # Endianness.qb contains newtype declarations that are needed more broadly,
    # so we include it first.
    "Endianness.qb",

    "TypeConversion.qb"

 
    "Math/NativeStubs.qb"
    "Math/Random.qb"

    


    "Paulis.qb",
    "ControlledOnBitString.qb"

    # # AmplitudeAmplification
    "AmplitudeAmplification/Utils.qb"
    "AmplitudeAmplification/Types.qb"
    "AmplitudeAmplification/AmplitudeAmplification.qb"
    
    #"AmplitudeAmplification/ExampleGrover.qb"
    "AmplitudeAmplification/ExampleAA.qb"   

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