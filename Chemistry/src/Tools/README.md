# Quantum Development Kit Chemistry Command-Line Tools

This project provides command-line tools for working with quantum chemistry data and the Quantum Development Kit chemistry library.

## Installing the Quantum Development Kit Chemistry Tools

The Quantum Development Kit chemistry command-line tools are distributed as a .NET Core SDK Global Tool, and thus can be installed using the `dotnet` command:

```
$ dotnet tool install --global Microsoft.Quantum.Chemistry.Tools
```

> ⚠ NOTE: If you want to use a pre-release version of the Microsoft.Quantum.Chemistry.Tools package, you'll need to [configure the Quantum Development Kit prerelease feed](https://github.com/microsoft/QuantumLibraries#optional-using-prerelease-versions) before proceeding.
> Once you have the prerelease feed configured, you can install the Microsoft.Quantum.Chemistry.Tools package by providing an explicit version:
> ```
> $ dotnet tool install --global Microsoft.Quantum.Chemistry.Tools --version "0.11.2005.1422-beta"
> ```

Once installed, you can run the Quantum Development Kit chemistry command-line tools with the `qdk-chem` command:

```
$ qdk-chem convert --from fcidump --to broombridge lih.fcidump --out lih.yaml
```

For help, use the `--help` option:

```
$ qdk-chem --help
$ qdk-chem convert --help
```

## Building from Source

Instead of installing the tools, you can also run `qdk-chem` directly from this project using the .NET Core SDK:

```
$ dotnet run -- convert --from fcidump --to broombridge lih.fcidump --out lih.yaml
```

Note that the `--` argument is used to separate arguments recognized by the .NET Core SDK from the arguments to the `qdk-chem` tool.

## Converting Chemistry Data with `qdk-chem`

The `qdk-chem convert` subcommand can be used to convert between the various chemistry data formats supported by the Quantum Development Kit.

> We strongly recommend the use of [Broombridge](https://docs.microsoft.com/quantum/libraries/chemistry/schema/broombridge) to represent quantum chemistry problems, and provide `qdk-chem convert` to make it easier to work with and prepare Broombridge-formatted data.

The source and destination formats can be specified by the `--from` and `--to` options, respectively.
By default, `qdk-chem` will write its output to the console; use the `--out` option to write directly to an output file.
If you want to provide the source file to be converted using stdin, pass `-` as the name of the file to read from:

```
$ echo lih.fcidump | qdk-chem convert --from fcidump --to broombridge - --out lih.yaml
```
