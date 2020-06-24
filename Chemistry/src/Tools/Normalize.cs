// Copyright (c) Microsoft Corporation. All rights reserved.
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

using Broombridge = Microsoft.Quantum.Chemistry.Broombridge;
using Microsoft.Quantum.Chemistry.Broombridge;

namespace Microsoft.Quantum.Chemistry.Tools
{

    public static class Normalize
    {
        public static Command CreateCommand() =>
            new Command("normalize")
            {
                new Argument<FileInfo>(
                    "path",
                    "Input data to be loaded, or - to load from stdin."
                ),
                new Option<SerializationFormat>(
                    "--format",
                    "Format to use in reading and writing problem description data."
                ),
                new Option<FileInfo?>(
                    "--out",
                    "Path to write output to. Data will be written to stdout by default."
                )
            }
            .WithDescription(
                "Normalizes the problem description data in a given problem description file" +
                "by loading and then saving the data in a given format."
            )
            .WithHandler<FileInfo, SerializationFormat, FileInfo?>(
                (path, format, @out) =>
                {
                    using var reader =
                        path.Name == "-"
                        ? System.Console.In
                        : File.OpenText(path.FullName);
                    using var writer =
                        @out == null
                        ? System.Console.Out
                        : new StreamWriter(File.OpenWrite(@out.FullName));
                    // Reuse logic from the convert command, but with the from
                    // and to format as the same.
                    Convert.ConvertProblemDescription(reader, format, writer, format);
                }
            );
    }   
 
}
