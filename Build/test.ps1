# Copyright (c) Microsoft Corporation. All rights reserved.
# Licensed under the MIT License.

$ErrorActionPreference = 'Stop'

& "$PSScriptRoot/set-env.ps1"
$all_ok = $True

function Test-One {
    Param($project)

    dotnet test (Join-Path $PSScriptRoot $project) `
        -c $Env:BUILD_CONFIGURATION `
        -v $Env:BUILD_VERBOSITY `
        --logger trx `
        /property:DefineConstants=$Env:ASSEMBLY_CONSTANTS `
        /property:Version=$Env:ASSEMBLY_VERSION

    return ($LastExitCode -eq 0) 
}

Write-Host "##[info]Testing Standard/tests/Standard.Tests.csproj"
$all_ok = (Test-One '../Standard/tests/Standard.Tests.csproj') -and $all_ok

Write-Host "##[info]Testing Chemistry/tests/ChemistryTests/QSharpTests.csproj"
$all_ok = (Test-One '../Chemistry/tests/ChemistryTests/QSharpTests.csproj') -and $all_ok

Write-Host "##[info]Testing Chemistry/tests/SystemTests/SystemTests.csproj"
$all_ok = (Test-One '../Chemistry/tests/SystemTests/SystemTests.csproj') -and $all_ok

Write-Host "##[info]Testing Chemistry/tests/DataModelTests/CSharpTests.csproj"
$all_ok = (Test-One '../Chemistry/tests/DataModelTests/CSharpTests.csproj') -and $all_ok

Write-Host "##[info]Testing Chemistry/tests/SerializationTests/SerializationTests.csproj"
$all_ok = (Test-One '../Chemistry/tests/SerializationTests/SerializationTests.csproj') -and $all_ok

Write-Host "##[info]Testing Chemistry/tests/JupyterTests/JupyterTests.csproj"
$all_ok = (Test-One '../Chemistry/tests/JupyterTests/JupyterTests.csproj') -and $all_ok

Write-Host "##[info]Testing Numerics/tests/NumericsTests.csproj"
$all_ok = (Test-One '../Numerics/tests/NumericsTests.csproj') -and $all_ok

if (-not $all_ok) {
    throw "At least one test failed execution. Check the logs."
}
