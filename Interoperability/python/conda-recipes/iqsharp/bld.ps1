# Constants
$DotnetVersion = "2.2.103";

# If the NUGET_VERSION environment variable is empty, default to PKG_VERSION.
if ($null -eq $Env:NUGET_VERSION) {
    $NugetVersion = $Env:PKG_VERSION;
} else {
    $NugetVersion = $Env:NUGET_VERSION;
}

# Convert conda architecture notation to dotnet notation.
if ($Env:ARCH -eq "32") {
    $DotnetArch = "x86";
} elseif ($Env:ARCH -eq "64") {
    $DotnetArch = "x64";
}

# Set environment variables to disable some .NET Core functionality that we don't want
# captured in the final package.
$Env:NUGET_XMLDOC_MODE = "skip"
$Env:DOTNET_SKIP_FIRST_TIME_EXPERIENCE = "1"

./dotnet-install.ps1 -Architecture $DotnetArch -Version $DotnetVersion -InstallDir $Env:LIBRARY_BIN -NoPath
# TODO: Change version of pkg below.
Get-Command jupyter | Write-Host

# Install IQ# itself.
dotnet tool install --tool-path $Env:LIBRARY_BIN --version "$NugetVersion"  Microsoft.Quantum.IQSharp --no-cache --add-source $Env:SRC_DIR

# Remove machine specific files from the .store folder for IQ#.
Remove-Item $Env:LIBRARY_BIN/.store/microsoft.quantum.iqsharp/*/project.assets.json
Remove-Item $Env:LIBRARY_BIN/.store/microsoft.quantum.iqsharp/*/restore.csproj.nuget.g.props

# Install IQ# into Jupyter
dotnet iqsharp install --prefix $Env:PREFIX
