#!/usr/bin/env pwsh
# Copyright (c) Microsoft Corporation.
# Licensed under the MIT License.

$quantumPackages = @(
    "Microsoft.Azure.Quantum.Client",
    "Microsoft.Azure.Quantum.DataStructures",
    "Microsoft.Azure.Quantum.Honeywell",
    "Microsoft.Azure.Quantum.IonQ",
    "Microsoft.Azure.Quantum.Machine",
    "Microsoft.Quantum.Chemistry",
    "Microsoft.Quantum.Chemistry.Jupyter",
    "Microsoft.Quantum.Compiler",
    "Microsoft.Quantum.CsharpGeneration",
    "Microsoft.Quantum.Decompositions.Honeywell",
    "Microsoft.Quantum.Decompositions.IonQ",
    "Microsoft.Quantum.Decompositions.OpenQASM",
    "Microsoft.Quantum.Decompositions.QCI",
    "Microsoft.Quantum.Development.Kit",
    "Microsoft.Quantum.IQSharp",
    "Microsoft.Quantum.IQSharp.Core",
    "Microsoft.Quantum.IQSharp.Jupyter",
    "Microsoft.Quantum.Katas",
    "Microsoft.Quantum.MachineLearning",
    "Microsoft.Quantum.Numerics",
    "Microsoft.Quantum.ProjectTemplates",
    "Microsoft.Quantum.Providers.Core",
    "Microsoft.Quantum.Providers.Datastructures",
    "Microsoft.Quantum.Providers.Honeywell",
    "Microsoft.Quantum.Providers.IonQ",
    "Microsoft.Quantum.Providers.QCI",
    "Microsoft.Quantum.QSharp.Core",
    "Microsoft.Quantum.Research",
    "Microsoft.Quantum.Research.Characterization",
    "Microsoft.Quantum.Research.Chemistry",
    "Microsoft.Quantum.Research.Simulation",
    "Microsoft.Quantum.Runtime.Core",
    "Microsoft.Quantum.Sdk",
    "Microsoft.Quantum.Simulators",
    "Microsoft.Quantum.Standard",
    "Microsoft.Quantum.Xunit"
);

function Update-QuantumProject() {
    <#
        .SYNOPSIS
            This cmdlet updates the Quantum SDK version as well as the references to Microsoft Quantum packages
            to a specified version of the QDK. If no version is specified, then the QDK version of the PowerShell
            Microsoft.Quantum.Utilities module is used. The dependencies of the project are restored, unless 
            parameter -NoRestore is specified.

        .PARAMETER Path
            Path to the VS project file to update.

        .PARAMETER Revert
            Restores a the target project to a previous version if available.
            Dependency restoration is not performed during a revert.

        .PARAMETER Version
            (Optional) Version of the QDK to update the project to. If not specified, the version of this PowerShell
            module will be used. To display the version of this modue, use cmdlet 'Get-QdkVersion'.

        .PARAMETER NoRestore
            (Optional) Skips performing the .NET restore step on the project after updating the references.

        .EXAMPLE
            Update-QuantumProject -Path QuantumFourierTransform.csproj

        .EXAMPLE
            Update-QuantumProject -Path QuantumFourierTransform.csproj -Version 0.11.2004.2825

        .EXAMPLE
            Update-QuantumProject -Path QuantumFourierTransform.csproj -Version 0.11.2004.2825 -NoRestore

        .EXAMPLE
            Update-QuantumProject -Path QuantumFourierTransform.csproj -Revert
            
    #>
    [CmdletBinding(SupportsShouldProcess=$true)]
    param(
        [string] $Version = $MyInvocation.MyCommand.Module.Version,

        [Parameter(Mandatory=$true, ValueFromPipeline=$true)]
        [System.IO.FileInfo] $Path,

        [Parameter(Mandatory=$false)]
        [switch] $Revert,

        [Parameter(Mandatory=$false)]
        [switch] $NoRestore
    );

    process {

        if (-Not (Test-Path $Path))
        {
            Write-Error "Project $Path does not exist.";
            return;
        }

        $prevVerDir = Join-Path -Path $Path.Directory -ChildPath ".prev";
        $pathToPrev = [System.IO.FileInfo](Join-Path -Path $prevVerDir -ChildPath $Path.Name);

        if ($Revert){
            Write-Verbose "Reverting $Path to previous version."

            if (Test-Path $pathToPrev) {
                if ($PSCmdlet.ShouldProcess($Path, "Revert")) {
                    Move-Item -Path $pathToPrev -Destination $Path -Force;
                }
            }
            else {
                Write-Output "A previous version of the project was not found. Skipping."
            }
        }
        else {
            Write-Verbose "Updating $Path to QDK version $Version."

            $csproj = [XML](Get-Content $Path);
            $projectWasUpdated = $false;

            # Update the Quantum SDK Version
            foreach ($item in (Select-Xml -Xml $csproj -XPath '//Project[@Sdk]')) {
                $sdk = [String]($item.node.Sdk);
                $sdk = $sdk -replace ("Microsoft.Quantum.Sdk`/`([^`"]*`)"),
                                     ("Microsoft.Quantum.Sdk`/" + $Version);

                if($sdk -ne $item.node.Sdk) {
                    $projectWasUpdated = $true;
                    $item.node.Sdk = $sdk;
                }
            }

            # Update the version of each of the Quantum packages
            foreach ($pack in $quantumPackages) {
                foreach ($item in (Select-Xml -Xml $csproj -XPath '//Project/ItemGroup/PackageReference')) {
                    if (($item.node.Include -ieq $pack) -and ($item.node.Version -ne $Version)) {
                        $item.node.Version = $Version;
                        $projectWasUpdated = $true;
                    }
                }
            }

            # Backup previous XML if a change was made and save updated file.
            if ($projectWasUpdated) {
                Write-Verbose "Saving previous configuration to $pathToPrev.";
                
                if ($PSCmdlet.ShouldProcess($pathToPrev, "Set-Content")) {
                    New-Item -ErrorAction Ignore -ItemType Directory -Path $prevVerDir;
                    Copy-Item -Path $Path -Destination $pathToPrev;
                }

                # Save the updated file.
                $csproj.Save($Path);
            }

            # Do a project restore of dependencies
            if (-Not $NoRestore) {
                dotnet restore $Path;
            }
        }
    }
}

function Get-QdkVersion() {
    <#
        .SYNOPSIS
            Displays the version of the Quantum Development Kit components.

        .DESCRIPTION
            Shows the version of the following components.
             - Microsoft.Quantum.Utilities
             - .NET Core SDK
             - Microsoft.Quantum.IQSharp
             - qsharp.py
             - VS Code Extension

        .EXAMPLE
            Get-QdkVersion
    #>

    [CmdletBinding()]
    param(
    );

    process {
        try{
            $dotnetVersion = (dotnet --version);
            $iqsharpVersion = [string]((dotnet iqsharp --version) -match "iqsharp" -replace "iqsharp:\s*","");
        }
        catch {
            # We'll just skip the values that couldn't be retrieved, and we will present empty strings
            # indicating that the version could not be retrieved. We only log the error in -Verbose mode.
            Write-Verbose $_;
        }

        try{
            $qsharpPythonVersion = (python -c "import qsharp; print(qsharp.__version__)");
        }
        catch {
            Write-Verbose $_;
        }

        @{
            "Microsoft.Quantum.Utilities" = $MyInvocation.MyCommand.Module.Version;
            ".NET Core SDK" = $dotnetVersion;
            "Microsoft.Quantum.IQSharp" = $iqsharpVersion;
            "qsharp.py" = $qsharpPythonVersion;
            "VS Code Extension" = (code --list-extensions --show-versions | Select-String quantum.quantum-devkit-vscode);
        }.GetEnumerator() | ForEach-Object { [pscustomobject]$_ }
    }
}
