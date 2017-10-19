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

function Find-QflatDependencies {
    [CmdletBinding()]
    param(
        [string]
        $Path
    )

    # Get a list of what we need to compile.
    # We do so by examining the *.csproj inside, looking for <None>
    # elements referencing *.qb files.
    # TODO: move to use <Compile> items with a DependentOn instead.
    # Import the csproj as an XML document.
    $csproj = [xml](Get-Content $Path)
    $csproj.Project.ItemGroup `
        | ForEach-Object { $_.None } `
        | Select-Object -ExpandProperty Include `
        | Where-Object { $_.EndsWith(".qb") } `
        | ForEach-Object { Join-Path $libDirectory $_ } `
        | Write-Output

}

function Find-QflatPrelude {
    [CmdletBinding()]
    param()

    Find-QflatCompiler `
        | Split-Path -Resolve `
        | ForEach-Object { Join-Path $_ ..\..\..\..\Library\standard.qb } `
        | Resolve-Path `
        | Write-Output
    
}
function Invoke-Qbc {
    param(
        [string] $Qbc,
        [string[]] $Sources
    )

    "Invoking qbc for files: " + $paths
    & $Qbc --input @sources;
    # $output = & $qbc --input @sources;

    # $stdout = $output | Where-Object { $_ -is [string]};
    # $stderr = $output | Where-Object { $_ -is [System.Management.Automation.ErrorRecord]};

    # Write-Verbose ($stdout -join "`n");
    if ($LASTEXITCODE -ne 0) {
        Write-Error "Qbc.exe failed:`n$($stderr -join "`n")";
        exit $LASTEXITCODE
    }

    #$output | Write-Output;
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