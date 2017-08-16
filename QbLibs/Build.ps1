##
# .SYNOPSIS
#     ./Build.ps1: This PowerShell script builds all QbLibs to generated C#
#     sources.
#
# .EXAMPLE
#     PS> $Env:QFLAT_PATH = "<path to Compiler.exe>"
#     PS> ./Build.ps1
##

# Set a default namespace.
# TODO: make this a parameter.
$Namespace = "Microsoft.Quantum.Canon"

function Out-GeneratedCSharp {
    [CmdletBinding()]
    param(
        [Parameter(Mandatory=$true, Position=1, ValueFromPipeline=$true)]
        [string[]] $FilePath
    )

    begin {
        # TODO: This will likely change to qbc.exe.
        $qbcBaseName = "Compiler.exe"

        # Find the compiler.
        if (Get-Command Compiler.exe -ErrorAction SilentlyContinue) {
            # It's on our $Env:PATH, so no need to worry.
            $qbc = $qbcBaseName
        } elseif ($Env:QFLAT_PATH) {
            # Try to find from the environment.
            $qbc = Join-Path $Env:QFLAT_PATH $qbcBaseName
        } else {
            # Make an educated guess that the user also has Solid checked out.
            # FIXME: This is a hack.
            $guessedPath = Join-Path "$PSScriptRoot/../../Solid/qflat/Compiler/bin/Debug/" $qbcBaseName;
            if (Test-Path $guessedPath) {
                $qbc = $guessedPath;
            } else {
                Write-Debug $guessedPath
                $msg = "Could not find $qbcBaseName."
                throw $msg;
            }
        }

        Write-Debug "Found compiler at $qbc."

        $targetDir = "$PSScriptRoot/Debug/$Namespace"

        if (!(Test-Path $targetDir)) {
            New-Item -ItemType Directory $targetDir
        }
    }

    process {
        $FilePath | ForEach-Object {
            $target = [IO.Path]::ChangeExtension((Join-Path $targetDir (Split-Path -Leaf $_)), "qb.generated.cs")
            Write-Debug "Generating code to $target."
            # TODO: make sure we put each generated C♯ file into the right subdir.
            & $qbc --operation generate --outputformat csharp --namespace $Namespace --input $_ --outputfile $target
        }
    }

    end {

    }
    
}

# Get a list of what we need to compile.
$libDirectory = Join-Path $PSScriptRoot $Namespace
$qflatSources = Get-ChildItem -Recurse (Join-Path $libDirectory "*.qb")
$qflatSources | Out-GeneratedCSharp