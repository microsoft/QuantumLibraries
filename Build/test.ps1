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

    if  ($LastExitCode -ne 0) {
        Write-Host "##vso[task.logissue type=error;]Failed to test $project."
        $script:all_ok = $False
    }
}

Write-Host "##[info]Testing Standard/tests/Standard.Tests.csproj"
Test-One '../Standard/tests/Standard.Tests.csproj'

Write-Host "##[info]Testing Chemistry/tests/ChemistryTests/QSharpTests.csproj"
Test-One '../Chemistry/tests/ChemistryTests/QSharpTests.csproj'

Write-Host "##[info]Testing Chemistry/tests/SystemTests/SystemTests.csproj"
Test-One '../Chemistry/tests/SystemTests/SystemTests.csproj'

Write-Host "##[info]Testing Chemistry/tests/DataModelTests/DataModelTests.csproj"
Test-One '../Chemistry/tests/DataModelTests/DataModelTests.csproj'

Write-Host "##[info]Testing Chemistry/tests/JupyterTests/JupyterTests.csproj"
Test-One '../Chemistry/tests/JupyterTests/JupyterTests.csproj'

Write-Host "##[info]Testing Numerics/tests/NumericsTests.csproj"
Test-One '../Numerics/tests/NumericsTests.csproj'

Write-Host "##[info]Testing QAOA/tests/QAOATests.csproj"
Test-One '../QAOA/tests/QAOATests.csproj'

if (-not $all_ok) {
    throw "At least one test failed execution. Check the logs."
}
