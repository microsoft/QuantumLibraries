function Find-QflatCompiler {
    [CmdletBinding()]
    param()

    $qbcBaseName = "qbc.exe"
    
    # Find the compiler.
    if (Get-Command $qbcBaseName -ErrorAction SilentlyContinue) {
        # It's on our $Env:PATH, so no need to worry.
        $qbc = $qbcBaseName
    } elseif ($Env:QFLAT_PATH) {
        # Try to find from the environment.
        $qbc = Join-Path $Env:QFLAT_PATH $qbcBaseName
    } else {
        # Make an educated guess that the user also has Solid checked out.
        # FIXME: This is a hack.
        $guessedPath = Join-Path "$PSScriptRoot/../../Solid/qflat/qbc/bin/Debug/" $qbcBaseName;
        if (Test-Path $guessedPath) {
            $qbc = $guessedPath;
        } else {
            Write-Debug $guessedPath
            $msg = "Could not find $qbcBaseName."
            throw $msg;
        }
    }

    Write-Debug "Found compiler at $qbc."

    $qbc | Write-Output;
}


function Invoke-Qbc {
    param(
        [string] $Qbc,
        [string[]] $Sources
    )

    $output = & $qbc --input @sources 2>&1;

    $stdout = $output | Where-Object { $_ -is [string]};
    $stderr = $output | Where-Object { $_ -is [System.Management.Automation.ErrorRecord]};

    Write-Verbose ($stdout -join "`n");
    if ($LASTEXITCODE -ne 0) {
        Write-Error "Qbc.exe failed:`n$($stderr -join "`n")";
    }
}

function ConvertFrom-Qflat {
    [CmdletBinding(
        DefaultParameterSetName="SingleFile"
    )]
    param(
        [Parameter(Mandatory=$true, Position=1, ValueFromPipeline=$true)]
        [string[]] $FilePath
    )

    begin {
        $qbc = Find-QflatCompiler;
        $paths = @()
    }

    process {
        $paths += @($FilePath)
    }

    end {
        Invoke-Qbc -Qbc $qbc -Sources $paths
    }
}