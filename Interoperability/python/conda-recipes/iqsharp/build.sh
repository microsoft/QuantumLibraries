# Constants
DOTNET_VERSION=2.2.103

# Convert conda architecture notation to dotnet notation.
if [ "$ARCH" == "32" ]; then
    DN_ARCH=x86
elif [ "$ARCH" == "64" ]; then
    DN_ARCH=x64
fi

# If the NUGET_VERSION environment variable is empty, default to PKG_VERSION.
if [ -z "$NUGET_VERSION" ]; then
    NUGET_VERSION=$PKG_VERSION
fi

# Set environment variables to disable some .NET Core functionality that we don't want
# captured in the final package.
export NUGET_XMLDOC_MODE=skip
export DOTNET_SKIP_FIRST_TIME_EXPERIENCE=1

# In order for .NET global tools to work correctly, we need to
# set an environment variable pointing to where they got installed.
# See https://github.com/dotnet/cli/issues/9114 for more details about DOTNET_ROOT.
export DOTNET_ROOT=$PREFIX/bin

echo "DOTNET_ROOT = $DOTNET_ROOT"

chmod +x ./dotnet-install.sh
./dotnet-install.sh --version $DOTNET_VERSION --architecture $DN_ARCH --install-dir $DOTNET_ROOT --no-path

# Install IQ# itself.
dotnet tool install --tool-path $PREFIX/bin --version $NUGET_VERSION Microsoft.Quantum.IQSharp --no-cache --add-source $SRC_DIR

# Remove machine specific files from the .store folder for IQ#.
rm -r $PREFIX/bin/.store/microsoft.quantum.iqsharp/*/project.assets.json
rm -r $PREFIX/bin/.store/microsoft.quantum.iqsharp/*/restore.csproj.nuget.g.props

# Install IQ# into Jupyter
dotnet iqsharp install --prefix $PREFIX

# Following the advice at https://conda-forge.org/docs/recipe.html, we'll need to make activate.sh and deactivate.sh
# commands in order to set the DOTNET_ROOT variable.
mkdir -p "${PREFIX}/etc/conda/activate.d"
printf "export _OLD_DOTNET_ROOT=\$DOTNET_ROOT\nexport DOTNET_ROOT=\$CONDA_PREFIX/bin" > "${PREFIX}/etc/conda/activate.d/${PKG_NAME}_activate.sh"
mkdir -p "${PREFIX}/etc/conda/deactivate.d"
printf "export DOTNET_ROOT=\$_OLD_DOTNET_ROOT" > "${PREFIX}/etc/conda/deactivate.d/${PKG_NAME}_deactivate.sh"
