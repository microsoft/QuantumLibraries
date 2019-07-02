# Copyright (c) Microsoft Corporation. All rights reserved.
# Licensed under the MIT License.

$ErrorActionPreference = 'Stop'

& "$PSScriptRoot/set-env.ps1"
$all_ok = $True

function Build-One {
    param(
        [string]$action,
        [string]$project
    );

    dotnet $action (Join-Path $PSScriptRoot $project) `
        -c $Env:BUILD_CONFIGURATION `
        -v $Env:BUILD_VERBOSITY `
        /property:DefineConstants=$Env:ASSEMBLY_CONSTANTS `
        /property:Version=$Env:ASSEMBLY_VERSION `
        /property:QsharpDocsOutDir=$Env:DOCS_OUTDIR

    return ($LastExitCode -eq 0)
}

Write-Host "##[info]Build Standard library"
$all_ok = (Build-One 'publish' '../Standard.sln') -and $all_ok

Write-Host "##[info]Build Chemistry library"
$all_ok = (Build-One 'publish' '../Chemistry.sln') -and $all_ok

Write-Host "##[info]Build Numerics library"
$all_ok = (Build-One 'publish' '../Numerics.sln') -and $all_ok

Write-Host "##[info]Build Standard library"
$all_ok = (Build-One 'publish' '../Magic.sln') -and $all_ok

if (-not $all_ok) {
    throw "At least one test failed execution. Check the logs."
}