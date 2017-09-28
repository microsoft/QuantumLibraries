##
# .SYNOPSIS
#     ./Build.ps1: This PowerShell script builds all QbLibs to generated C#
#     sources.
#
# .EXAMPLE
#     PS> $Env:QFLAT_PATH = "<path to Compiler.exe>"
#     PS> ./Build.ps1
##
param(
    [switch]
    $Batch = $false,

    [switch]
    $ParseOnly = $false
)

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

# Get a list of what we need to compile.
# Unlike before, we manually specify since the order matters.
$libDirectory = Join-Path $PSScriptRoot "Microsoft.Quantum.Canon"
$qflatSources = @(
    # Provide stubs for primitive operations.
    "Stubs.qb",

    # Endianness.qb contains newtype declarations that are needed more broadly,
    # so we include it first.
    "Endianness.qb",

    # Similarly with OracleTypes.qb, save for that it depends on OperationPow.qb.
    "OperationPow.qb",
    "OracleTypes.qb",

    "ApplyToEach.qb",
    "ApplyToRange.qb",
    "Arithmetic.qb",
    "IterativePhaseEstimation.qb",
    # # "QFT.qb", # QFT commented out in lieu of merging in martinro/ branch.
    "QuantumPhaseEstimation.qb",
    "ShiftOp.qb",
    "With.qb",

    "Paulis.qb",
    
    # # QECC
    "Qecc/Types.qb",
    "Qecc/Utils.qb",
    "Qecc/BitFlipCode.qb"
) | ForEach-Object {
    Join-Path $libDirectory $_
}

if ($Batch) {
    $qflatSources | ConvertFrom-Qflat -ParseOnly:$ParseOnly -Batch -Target Canon.g.cs
} else {
    $qflatSources | ConvertFrom-Qflat -ParseOnly:$ParseOnly
}