#!/bin/bash 

ver=$1
pkgs=$2

# It is recommended that you clean up the environment first
# by running this command:
# git submodule deinit --all -f
git submodule update --init

if [ "$(uname)" == "Darwin" ]; then
  backup="''"
else
  backup=""
fi

: ${ver:="$NUGET_VERSION"}
: ${pkgs:="Microsoft.Quantum.Development.Kit;Microsoft.Quantum.IQSharp.Core;Microsoft.Quantum.Simulators;Microsoft.Quantum.Compiler;Microsoft.Quantum.Canon;Microsoft.Quantum.Xunit;Microsoft.Quantum.Chemistry;Microsoft.Quantum.Research"}


for pkg in `echo $pkgs | tr ";" "\n"`; do 
  echo Will update package $pkg with version $ver...

  grep --include=\packages.config -lri -e "package *id=\"$pkg\" *version=" * | xargs sed -i $backup "s/package *id=\"$pkg\" *version=\"\([^\"]*\)\"/package id=\"$pkg\" version=\"$ver\"/i"
  grep --include=\*proj -lri -e "PackageReference *Include=\"$pkg\" *Version=" * | xargs sed -i $backup "s/PackageReference *Include=\"$pkg\" *Version=\"\([^\"]*\)\"/PackageReference Include=\"$pkg\" Version=\"$ver\"/i"
done 

echo done!
echo

# Commit changes to a release branch:
git status
git --no-pager diff

git submodule foreach git status
git submodule foreach git --no-pager diff
git submodule foreach git checkout -b release/v$ver

echo 
echo ------------------------------------------------------------------
echo - Changes not pushed to origin; to push changes run this command: -
echo ------------------------------------------------------------------
echo git submodule foreach git commit --all -m \"Updating projects to $ver\" --allow-empty
echo git submodule foreach git push --set-upstream origin to release/v$ver
