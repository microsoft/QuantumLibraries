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
        /property:QsharpDocsOutputPath=$Env:DOCS_OUTDIR

    if  ($LastExitCode -ne 0) {
        Write-Host "##vso[task.logissue type=error;]Failed to build $project."
        $script:all_ok = $False
    }
}

function Build-PwshDocs {
    param(
        [string] $DocsPath,
        [string] $ModulePath
    );

    # Check that platyPS is available; if not, we can't
    # build docs.
    if (!(Get-Module -ListAvailable platyPS)) {
        Write-Host "##vso[task.logissue type=warning;]platyPS was not available, cannot build docs for PowerShell modules.";
    } else {
        try {
            Import-Module platyPS
            $platyPS = Get-Module platyPS;
            Write-Host "##[info]Using platyPS version $($platyPS.Version) to build docs.";
            New-ExternalHelp -Force $DocsPath -OutputPath (Join-Path $ModulePath "en-US") -Verbose;
        } catch {
            Write-Verbose "$_"
            $script:all_ok = $false;
        }
    }
}

Write-Host "##[info]Build Standard library"
Build-One 'publish' '../Standard.sln'

Write-Host "##[info]Build Chemistry library"
Build-One 'publish' '../Chemistry.sln'

Write-Host "##[info]Build Numerics library"
Build-One 'publish' '../Numerics.sln'

Write-Host "##[info]Build QML library"
Build-One 'publish' '../MachineLearning.sln'

Write-Host "##[info]Build Jupyter magic library"
Build-One 'publish' '../Magic.sln'

Write-Host "##[info]Building PowerShell docs..."
Build-PwshDocs `
    -ModulePath (Join-Path $PSScriptRoot "Utilities/src/PowerShell/Microsoft.Quantum.Utilities") `
    -DocsPath (Join-Path $PSScriptRoot "Utilities/docs")

if (-not $all_ok) {
    throw "At least one test failed execution. Check the logs."
}
