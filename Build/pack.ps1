# Copyright (c) Microsoft Corporation. All rights reserved.
# Licensed under the MIT License.

$ErrorActionPreference = 'Stop'

& "$PSScriptRoot/set-env.ps1"
$all_ok = $True

function Pack-One() {
    Param($project)

    dotnet pack (Join-Path $PSScriptRoot $project) `
        --no-build `
        -c $Env:BUILD_CONFIGURATION `
        -v $Env:BUILD_VERBOSITY `
        -o $Env:NUGET_OUTDIR `
        /property:PackageVersion=$Env:NUGET_VERSION

    if  ($LastExitCode -ne 0) {
        Write-Host "##vso[task.logissue type=error;]Failed to pack $project."
        $script:all_ok = $False
    }
}

function Pack-Wheel() {
    param(
        [string] $Path
    );

    $result = 0

    Push-Location (Join-Path $PSScriptRoot $Path)
        python setup.py bdist_wheel sdist --formats=gztar

        if  ($LastExitCode -ne 0) {
            Write-Host "##vso[task.logissue type=error;]Failed to build $Path."
            $script:all_ok = $False
        } else {
            Copy-Item "dist/*.whl" $Env:PYTHON_OUTDIR
            Copy-Item "dist/*.tar.gz" $Env:PYTHON_OUTDIR
        }
    Pop-Location
}


Write-Host "##[info]Pack Standard library"
Pack-One '../Standard/src/Standard.csproj'

Write-Host "##[info]Pack Standard visualization library"
Pack-One '../Visualization/src/Visualization.csproj'

Write-Host "##[info]Pack Chemistry library"
Pack-One '../Chemistry/src/Runtime/Runtime.csproj'
Pack-One '../Chemistry/src/DataModel/DataModel.csproj'
Pack-One '../Chemistry/src/Metapackage/Metapackage.csproj'

Write-Host "##[info]Pack QML library"
Pack-One '../MachineLearning/src/MachineLearning.csproj'

Write-Host "##[info]Pack Numerics library"
Pack-One '../Numerics/src/Numerics.csproj'

Write-Host "##[info]Pack chemistry magics library"
Pack-One '../Chemistry/src/Jupyter/Jupyter.csproj'

Write-Host "##[info]Pack chemistry tool"
Pack-One '../Chemistry/src/Tools/Tools.csproj'

if ($Env:ENABLE_PYTHON -eq "false") {
    Write-Host "##vso[task.logissue type=warning;]Skipping Creating Python packages. Env:ENABLE_PYTHON was set to 'false'."
} else {
    Write-Host "##[info]Packing Python wheel..."
    python --version
    Pack-Wheel '../Python/qsharp-chemistry'
    Pack-Wheel '../Python/qsharp'
}

if (-not $all_ok) {
    throw "At least one test failed execution. Check the logs."
}
