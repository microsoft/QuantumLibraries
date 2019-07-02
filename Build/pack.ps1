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

    return ($LastExitCode -eq 0)
}

Write-Host "##[info]Pack Standard library"
$all_ok = (Pack-One '../Standard/src/Standard.csproj') -and $all_ok

Write-Host "##[info]Pack Chemistry library"
$all_ok = (Pack-One '../Chemistry/src/DataModel/DataModel.csproj') -and $all_ok

Write-Host "##[info]Pack Numerics library"
$all_ok = (Pack-One '../Numerics/src/Numerics.csproj') -and $all_ok

Write-Host "##[info]Pack Standard library"
$all_ok = (Pack-One '../Chemistry/src/Jupyter/Jupyter.csproj') -and $all_ok

if (-not $all_ok) {
    throw "At least one test failed execution. Check the logs."
}