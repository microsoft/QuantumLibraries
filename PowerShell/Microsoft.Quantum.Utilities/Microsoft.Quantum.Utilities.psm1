function Update-QuantumProject() {
    [CmdletBinding()]
    param(
        [string] $Version = $MyInvocation.MyCommand.Module.Version,

        [Parameter(Mandatory=$true, ValueFromPipeline=$true)]
        [System.IO.FileInfo] $Path
    );

    begin {
        Write-Verbose "Updating project $Path to QDK $Version...";
    }

    process {
        Write-Output "OK";
    }
}

function Get-QdkVersion() {
    [CmdletBinding()]
    param(
    );

    process {
        Write-Output $MyInvocation.MyCommand.Module.Version;
    }
}