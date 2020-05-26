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

function Test-PowerShellModules {
    param(
        [string]
        $Path
    );

    if (!(Get-Module -ListAvailable Pester)) {
        Write-Host "##vso[task.logissue type=warning;]Pester not available, cannot run PowerShell tests.";
    } else {
        Import-Module Pester
        $results = Invoke-Pester $Path -PassThru;
        $Script:all_ok = $Script:all_ok -and ($results.FailedCount -eq 0);
    }
}

Write-Host "##[info]Testing Standard/tests/Standard.Tests.csproj"
Test-One '../Standard/tests/Standard.Tests.csproj'

Write-Host "##[info]Testing Chemistry/tests/ChemistryTests/QSharpTests.csproj"
Test-One '../Chemistry/tests/ChemistryTests/QSharpTests.csproj'

Write-Host "##[info]Testing Chemistry/tests/SystemTests/SystemTests.csproj"
Test-One '../Chemistry/tests/SystemTests/SystemTests.csproj'

Write-Host "##[info]Testing Chemistry/tests/DataModelTests/CSharpTests.csproj"
Test-One '../Chemistry/tests/DataModelTests/CSharpTests.csproj'

Write-Host "##[info]Testing Chemistry/tests/SerializationTests/SerializationTests.csproj"
Test-One '../Chemistry/tests/SerializationTests/SerializationTests.csproj'

Write-Host "##[info]Testing Chemistry/tests/JupyterTests/JupyterTests.csproj"
Test-One '../Chemistry/tests/JupyterTests/JupyterTests.csproj'

Write-Host "##[info]Testing Numerics/tests/NumericsTests.csproj"
Test-One '../Numerics/tests/NumericsTests.csproj'

if (-not $all_ok) {
    throw "At least one test failed execution. Check the logs."
}
