# Microsoft Quantum Development Kit Libraries #

Welcome to the Microsoft Quantum Development Kit!

This repository contains open-source libraries for the [Quantum Development Kit](https://docs.microsoft.com/azure/quantum):

- **[Docs/](./Docs)**: Additional documentation for developing on the libraries. Please see [QDK online documentation](https://docs.microsoft.com/azure/quantum/) for online documentation.
- **[Standard/](./Standard)**: Q# sources used to implement [the Q# standard libraries](https://docs.microsoft.com/azure/quantum/user-guide/libraries/standard).
- **[Chemistry/](./Chemistry)**: Q# and C# sources used to implement a library for [quantum chemistry](https://docs.microsoft.com/azure/quantum/user-guide/libraries/chemistry) and Hamiltonian simulation.
- **[Numerics/](./Numerics)**: Q# sources used to implement the [quantum numerics library](https://docs.microsoft.com/azure/quantum/user-guide/libraries/numerics/numerics).
- **[LICENSE](./LICENSE.txt)**: Terms of use and license details for the Quantum Development Kit libraries.

## New to Quantum? ##

See the [introduction to quantum computing](https://docs.microsoft.com/azure/quantum/concepts-overview/) provided with the Quantum Development Kit.

## Getting Started ##

The libraries provided in this repository are built using [.NET Core](https://docs.microsoft.com/en-us/dotnet/core/) and the
[Quantum Development Kit](https://docs.microsoft.com/azure/quantum).
Please see the [installation guide](https://docs.microsoft.com/azure/quantum/install-overview-qdk) for how to get up and running.

You may also visit our [Quantum](https://github.com/Microsoft/Quantum) repository, which offers a wide variety
of samples on how to use these libraries to write quantum based programs.

## Build Status ##

| branch | status    |
|--------|-----------|
| main | [![Build Status](https://dev.azure.com/ms-quantum-public/Microsoft%20Quantum%20(public)/_apis/build/status/Microsoft.QuantumLibraries?branchName=main)](https://dev.azure.com/ms-quantum-public/Microsoft%20Quantum%20(public)/_build/latest?definitionId=1&branchName=main) |

## Feedback ##

If you have feedback about the content in this repository, please let us know by filing a [new issue](https://github.com/microsoft/quantumlibraries/issues/new/choose)!
If you have feedback about some other part of the Microsoft Quantum Development Kit, please see the [contribution guide](https://docs.microsoft.com/azure/quantum/contributing-overview/) for more information.

## Contributing ##

This project welcomes contributions and suggestions.  Most contributions require you to agree to a
Contributor License Agreement (CLA) declaring that you have the right to, and actually do, grant us
the rights to use your contribution. For details, visit https://cla.microsoft.com.

When you submit a pull request, a CLA-bot will automatically determine whether you need to provide
a CLA and decorate the PR appropriately (e.g., label, comment). Simply follow the instructions
provided by the bot. You will only need to do this once across all repos using our CLA.

This project has adopted the [Microsoft Open Source Code of Conduct](https://opensource.microsoft.com/codeofconduct/).
For more information see the [Code of Conduct FAQ](https://opensource.microsoft.com/codeofconduct/faq/) or
contact [opencode@microsoft.com](mailto:opencode@microsoft.com) with any additional questions or comments.

## [Optional] Using Prerelease Versions ##

If you're interested in helping test the Quantum Development Kit libraries, or if you want to try out new features before they are released, you can add the [Quantum Development Kit prerelease feed](https://dev.azure.com/ms-quantum-public/Microsoft%20Quantum%20(public)/_packaging?_a=feed&feed=alpha) to your .NET Core SDK configuration.
Packages on the prerelease feed are marked with `-alpha` in their version number, so that projects built using released versions of Quantum Development Kit libraries will not be affected.
Note that the prerelease feed is used automatically when building libraries in this repository.

To use the prerelease feed, edit your `NuGet.Config` file to include the prerelease feed URL (`https://pkgs.dev.azure.com/ms-quantum-public/Microsoft Quantum (public)/_packaging/alpha/nuget/v3/index.json`) as a package source.
The location of this file varies depending on your operating system:

| OS | NuGet config file location |
|----|----------------------------|
| Windows | `$Env:APPDATA/Roaming/NuGet/NuGet.Config` |
| macOS / Linux | `~/.config/NuGet/NuGet.Config` or `~/.nuget/NuGet/NuGet.Config` |

Note that this file may not already exist, depending on your configuration.

For example, the following `NuGet.Config` file includes both the main NuGet package feed, and the Quantum Development Kit prerelease feed:
## Site Policy on GitHub
```xml
<?xml version="1.0" encoding="utf-8"?>
<configuration>
  <packageSources>
    <add key="nuget.org" value="https://api.nuget.org/v3/index.json" protocolVersion="3" />
    <add key="qdk-alpha" value="https://pkgs.dev.azure.com/ms-quantum-public/Microsoft Quantum (public)/_packaging/alpha/nuget/v3/index.json" protocolVersion="3" />
  </packageSources>
</configuration>
```
