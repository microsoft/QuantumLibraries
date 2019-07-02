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

    $script:all_ok = ($LastExitCode -eq 0) -and $script:all_ok
}

Write-Host "##[info]Build Standard library"
Build-One 'publish' '../Standard.sln'

Write-Host "##[info]Build Chemistry library"
Build-One 'publish' '../Chemistry.sln'

Write-Host "##[info]Build Numerics library"
Build-One 'publish' '../Numerics.sln'

Write-Host "##[info]Build Standard library"
Build-One 'publish' '../Magic.sln'

if (-not $all_ok) {
    throw "At least one test failed execution. Check the logs."
}