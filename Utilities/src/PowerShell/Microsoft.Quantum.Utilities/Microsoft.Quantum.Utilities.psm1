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
    "Microsoft.Quantum.IQSharp.Core",
    "Microsoft.Quantum.IQSharp.Jupyter",
    "Microsoft.Quantum.Numerics",
    "Microsoft.Quantum.QSharp.Core",
    "Microsoft.Quantum.Research",
    "Microsoft.Quantum.Runtime.Core",
    "Microsoft.Quantum.Simulators",
    "Microsoft.Quantum.Providers.Core",
    "Microsoft.Quantum.Providers.Datastructures",
    "Microsoft.Quantum.Providers.Honeywell",
    "Microsoft.Quantum.Providers.IonQ",
    "Microsoft.Quantum.Providers.QCI",
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
                    New-Item -ItemType Directory -Path $prevVerDir;
                    Move-Item -Path $pathToPrev -Destination $Path -Force;
                }
            }
        }
        else {
            Write-Verbose "Updating $Path to QDK version $Version."

            $csproj = Get-Content $Path;

            # Update the Quantum SDK Version
            $upcsproj = $csproj -replace ("Sdk=`"Microsoft.Quantum.Sdk`/`([^`"]*`)`""),
                                         ("Sdk=`"Microsoft.Quantum.Sdk`/" + $Version + "`"");

            # Update the version of each of the Quantum packages
            foreach ($pack in $quantumPackages) {
                $upcsproj  = $upcsproj -replace ("PackageReference\s*Include=`"" + $pack + "`"\s*Version=`"`([^`"]*`)`"") , 
                                                ("PackageReference Include=`"" + $pack + "`" Version=`"" + $Version + "`"");
            }

            Set-Content -Path $Path -Value $upcsproj;

            # Do a project restore of dependencies
            if (-Not $NoRestore) {
                dotnet restore $Path;
            }

            if (Compare-Object -ReferenceObject $csproj -DifferenceObject $upcsproj) {
                Write-Verbose "Saving previous configuration to $pathToPrev.";
                
                if ($PSCmdlet.ShouldProcess($pathToPrev, "Set-Content")) {
                    New-Item -ErrorAction Ignore -ItemType Directory -Path $prevVerDir;
                    Set-Content -Path $pathToPrev -Value $csproj;
                }
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
