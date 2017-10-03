function Split-Errors {
    param(
        [string[]] $Output,

        [string[]] $SplitAt = @("The following", "FAILURE: parsing", "resulted in"),

        [int] $Offset = 1
    )

    $LineNums = ($Output | Select-String $SplitAt -Context 0 | Select-Object -ExpandProperty LineNumber);
    $LineNum = $LineNums[$LineNums.Length - 1];
    $LineNum | write-host;

    # We need an extra -1 in calculating $Offset to deal with that Select-String uses 1-based indexing.
    ($Output[($LineNum + $Offset - 1)..$Output.Length] | % { "    " + $_ }) -join "`n"
}

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
        [string] $qbc,
        [string[]] $sources,
        [string] $target,
        [switch] $ParseOnly = $false
    )

    if ($ParseOnly) {
        $output = & $qbc --operation parse --outputformat csharp --input @sources --outputfile $target 2>&1;
    } else {
        $output = & $qbc --operation generate --outputformat csharp --input @sources --outputfile $target 2>&1;
    }

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
        [Parameter(Mandatory=$true, Position=1, ValueFromPipeline=$true, ParameterSetName="SingleFile")]
        [Parameter(Mandatory=$true, Position=1, ValueFromPipeline=$true, ParameterSetName="Batch")]
        [string[]] $FilePath,

        [Parameter(ParameterSetName="SingleFile")]
        [Parameter(ParameterSetName="Batch")]
        [switch] $ParseOnly = $false,
                
        [Parameter(ParameterSetName="Batch")]
        [switch] $Batch = $false,

        [Parameter(Position=2, ParameterSetName="Batch")]
        [string] $Target
    )

    begin {
        $qbc = Find-QflatCompiler;
        $batchMode = $PSCmdlet.ParameterSetName -eq "Batch";

        if ($batchMode) {
            $batchFiles = @();
        }
    }

    process {
        if ($batchMode) {
            $batchFiles += $FilePath;
        } else {
            $FilePath | ForEach-Object {
                $target = [IO.Path]::ChangeExtension($_, "g.cs")
                Write-Debug "Generating code to $target."
                Invoke-Qbc -qbc $qbc -sources @($_) -target $target -ParseOnly:$ParseOnly
            }
        }
    }

    end {
        if ($batchMode) {
            # TODO: deduplicate code here.
            Write-Debug "Generating code to $Target."
            Invoke-Qbc -qbc $qbc -sources $batchFiles -target $Target -ParseOnly:$ParseOnly
        }
    }
    
}