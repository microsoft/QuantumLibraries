// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

#nullable enable

using Microsoft.Extensions.Configuration;
using Microsoft.Extensions.Logging;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using System.CommandLine;
using System.CommandLine.Invocation;
using System.IO;
using System;

namespace Microsoft.Quantum.Chemistry.Tools
{
    
    public class Program
    {

        public static Command CreateRootCommand() =>
            new RootCommand
            {
                Convert.CreateCommand(),
                Normalize.CreateCommand(),
                ExportJW.CreateCommand()
            }
            .WithDescription("Tools for working with quantum chemistry data.");

        public static async Task<int> Main(string[] args) =>
            await CreateRootCommand().InvokeAsync(args);

    }
}
