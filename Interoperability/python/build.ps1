param(
    [string] $Config = $env:BuildConfiguration,
    [string] $PythonVersion = "3.5"
);


# Check conda installation, fail if missing.
if (-not (Get-Command conda -ErrorAction SilentlyContinue)) {
    Write-Error "conda (https://conda.io/docs/index.html) is required to build the qsharp package."
    exit 1;
}

# Check for conda-build, installed if missing.
if (-not (Get-Command conda-build -ErrorAction SilentlyContinue)) {
    Write-Host "Installing conda-build."
    conda install -y conda-build
}

# Make sure the ASSEMBLY_VERSION env is set:
if ("$env:ASSEMBLY_VERSION" -eq "") {
    Write-Warning "Missing ASSEMBLY_VERSION, defaulting to '0.0.0.0'"
    $env:ASSEMBLY_VERSION = '0.0.0.0'
}

# Make sure the Config env is set:
if ("$Config" -eq "") {
    Write-Warning "Missing BuildConfiguration, default to Release"
    $Config = 'Release'
}

# Report versions.
conda --version | Write-Host
conda-build --version | Write-Host
"Version $env:ASSEMBLY_VERSION" | Write-Host
"Config $Config" | Write-Host

# Build Canon and copy its dll to the package folder:
if (-not (Get-Item -Path ..\..\Canon\src\bin\$Config\netstandard2.0\Microsoft.Quantum.Canon.dll -ErrorAction SilentlyContinue)) {
    Write-Warning "Building Canon"
    dotnet build ..\..\Canon\src -c $Config /p:Version=$env:ASSEMBLY_VERSION
}
Write-Host "Copying Canon"
Copy-Item  -Path ..\..\Canon\src\bin\$Config\netstandard2.0\Microsoft.Quantum.*.dll -Destination qsharp

# Pre-build the package, so build's Nugets are picked up correctly.
Write-Host "Setting up package"
python setup.py build

# Use conda-build to build the package itself.
Write-Host "Calling conda-build"
conda-build . -c pythonnet -c conda-forge --python $PythonVersion --output-folder bin\$Config\package

