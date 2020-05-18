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
            This cmdlet upates the Quantum SDK version as well as the references to Microsoft Quantum packages
            to a specified version of the QDK. If no version is specified, then the QDK version of the PowerShell
            Microsoft.Quantum.Utilities module is used.

        .PARAMETER Path
            Path to the VS project file to update.

        .PARAMETER Revert
            Restores a the target project to a previous version if available. If used, parameter Version is ignored.

        .PARAMETER Version
            (Optional) Version of the QDK to update the project to. If not specified, the version of this PowerShell
            module will be used. To display the version of this modue, use cmdlet 'Get-QdkVersion'.

        .EXAMPLE
            Update-QuantumProject -Path QuantumFourierTransform.csproj

        .EXAMPLE
            Update-QuantumProject -Path QuantumFourierTransform.csproj -Version 0.11.2004.2825

        .EXAMPLE
            Update-QuantumProject -Path QuantumFourierTransform.csproj -Revert
            
    #>
    [CmdletBinding()]
    param(
        [string] $Version = $MyInvocation.MyCommand.Module.Version,

        [Parameter(Mandatory=$true, ValueFromPipeline=$true)]
        [System.IO.FileInfo] $Path,

        [Parameter(Mandatory=$false)]
        [switch] $Revert
    );

    process {

        if (-Not (Test-Path $Path))
        {
            Write-Error "Project $Path does not exist.";
            return;
        }

        $pathToPrev = [System.IO.FileInfo] ([String]$Path + ".prev_ver");

        if ($Revert){
            Write-Verbose "Reverting $Path to previous version."

            if (Test-Path $pathToPrev) {
                Move-Item -Path $pathToPrev -Destination $Path -Force;
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

            if (Compare-Object -ReferenceObject $csproj -DifferenceObject $upcsproj) {
                
                Set-Content -Path $pathToPrev -Value $csproj;
                Write-Verbose "Saving previous configuration to $pathToPrev.";
            }
        }
    }
}

function Get-QdkVersion() {
    <#
        .SYNOPSIS
            Displays the version of the Quantum Development Kit corresponding to this module.

        .EXAMPLE
            Get-QdkVersion
            
    #>

    [CmdletBinding()]
    param(
    );

    process {
        Write-Output $MyInvocation.MyCommand.Module.Version;
    }
}