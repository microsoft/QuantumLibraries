# Canon #

This folder contains draft operations for a Q♭ standard library, each of which is itself written in Q♭.
The operations implemented as a part of Canon should ideally be useful for combining user-defined operations into a larger quantum algorithm.

## Building Canon ##

Canon is provided with a PowerShell script, ``.\Build.ps1``, that locates and runs the Q♭ compiler to generate C♯ code from the Canon source.
The build script first checks if ``Compiler.exe`` is on ``$Env:PATH``, then checks the environment variable ``$Env:QFLAT_PATH``, and finally guesses that the Solid repository is checked out as a neighbor to the Libraries repo.

In any case, to build Canon, simply run the build script:
```powershell
PS> Unblock-File .\Build.ps1 # May be required, depending on Get-ExecutionPolicy.
PS> .\Build.ps1
```

**TODO:** port to MSBuild.

