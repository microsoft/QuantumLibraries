#!/usr/bin/env pwsh
#Requires -PSEdition Core

& "$PSScriptRoot/set-env.ps1"

@{
    Packages = @(
        "Microsoft.Quantum.Standard",
        "Microsoft.Quantum.Chemistry",
        "Microsoft.Quantum.Numerics",
        "Microsoft.Quantum.MachineLearning"
    );
    Assemblies = @(
        ".\Standard\src\bin\$Env:BUILD_CONFIGURATION\netstandard2.1\Microsoft.Quantum.Standard.dll",
        ".\Numerics\src\bin\$Env:BUILD_CONFIGURATION\netstandard2.1\Microsoft.Quantum.Numerics.dll",
        ".\MachineLearning\src\bin\$Env:BUILD_CONFIGURATION\netstandard2.1\Microsoft.Quantum.MachineLearning.dll",
        ".\Chemistry\src\DataModel\bin\$Env:BUILD_CONFIGURATION\netstandard2.1\Microsoft.Quantum.Chemistry.DataModel.dll",
        ".\Chemistry\src\Runtime\bin\$Env:BUILD_CONFIGURATION\netstandard2.1\Microsoft.Quantum.Chemistry.Runtime.dll"
    ) | ForEach-Object { Get-Item (Join-Path $PSScriptRoot ".." $_) };
} | Write-Output;
