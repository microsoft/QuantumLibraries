// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.
#nullable enable

using System.Collections.Generic;
using System.Linq;
using System.CommandLine;
using System.IO;
using System;

using Microsoft.Quantum.Chemistry.Broombridge;
using Microsoft.Quantum.Chemistry.OrbitalIntegrals;
using Microsoft.Quantum.Chemistry.Fermion;
using Microsoft.Quantum.Chemistry.QSharpFormat;
using System.Text.Json;
using System.Text.Json.Serialization;
using Microsoft.Quantum.Simulation.Core;

namespace Microsoft.Quantum.Chemistry.Tools
{

    public static class ExportJW
    {
        public static Command CreateCommand() =>
            new Command("export-jw")
            {
                new Argument<FileInfo>(
                    "path",
                    "Input data to be loaded, or - to load from stdin."
                ),
                new Option<SerializationFormat>(
                    "--format",
                    "Format to use in loading problem description data."
                ),
                new Option<FileInfo?>(
                    "--out",
                    "Path to write output to. Data will be written to stdout by default."
                ),
                new Option<bool>(
                    "--flatten",
                    "If true, flattens resulting JSON (often easier for use in " +
                    "native code)."
                )
            }
            .WithDescription(
                "Exports a JSON representation of the Jordan–Wigner transformation of " +
                "the fermionic Hamiltonian for a particular electronic structure problem."
            )
            .WithHandler<FileInfo, SerializationFormat, FileInfo?, bool>(
                (path, from, @out, flatten) =>
                {
                    using var reader =
                        path.Name == "-"
                        ? System.Console.In
                        : File.OpenText(path.FullName);
                    using var writer =
                        @out == null
                        ? System.Console.Out
                        : new StreamWriter(File.OpenWrite(@out.FullName));
                    ExportJwData(reader, from, writer, flatten: flatten);
                }
            );

        public static void ExportJwData(
            TextReader reader, SerializationFormat from,
            TextWriter writer,
            IndexConvention indexConvention = IndexConvention.UpDown,
            bool flatten = true
        )
        {
            var data = Load(reader, from).ToList();
            if (data.Count != 1)
            {
                System.Console.Error.WriteLine($"Expected a single problem description, but got a list of {data.Count}.");
            }
            var problem = data.Single();

            var fermionHamiltonian = problem
                .OrbitalIntegralHamiltonian
                .ToFermionHamiltonian(indexConvention);
            var jwHamiltonian = fermionHamiltonian
                .ToPauliHamiltonian(Paulis.QubitEncoding.JordanWigner)
                .ToQSharpFormat();
            var wavefunction = (
                    (problem.InitialStates?.Count ?? 0) == 0
                    ? fermionHamiltonian.CreateHartreeFockState(problem.NElectrons)
                    : problem
                        .InitialStates
                        .First()
                        .Value
                        .ToIndexing(indexConvention)
                )
                .ToQSharpFormat();

            // var encoded = JsonSerializer.Serialize(
            //     QSharpFormat.Convert.ToQSharpFormat(jwHamiltonian, wavefunction),
            //     options
            // );
            var qsData = QSharpFormat.Convert.ToQSharpFormat(jwHamiltonian, wavefunction);
            var encoded = flatten
                          ? JsonSerializer.Serialize(Flatten(qsData))
                          : JsonSerializer.Serialize(qsData);

            writer.Write(encoded);
            writer.Close();
        }

        // TODO: Move into common class.
        internal static IEnumerable<ElectronicStructureProblem> Load(TextReader reader, SerializationFormat from) =>
            (from switch
            {
                SerializationFormat.Broombridge =>
                    BroombridgeSerializer.Deserialize(reader),
                SerializationFormat.LiQuiD => 
                    LiQuiDSerializer.Deserialize(reader),
                SerializationFormat.FciDump =>
                    FciDumpSerializer.Deserialize(reader),
                _ => throw new ArgumentException($"Invalid format {from}.")
            })
            .ToList();

        internal static object[] Flatten(JordanWigner.JordanWignerEncodingData data) =>
            new object[]
            {
                data.Item1,
                Flatten(data.Item2),
                new object[]
                {
                    data.Item3.Item1,
                    data.Item3.Item2.Select(s => Flatten(s)).ToArray()
                },
                data.Item4
            };

        internal static object[] Flatten(JordanWigner.JWOptimizedHTerms data) =>
            new object[]
            {
                data.Item1.Select(s => Flatten(s)).ToArray(),
                data.Item2.Select(s => Flatten(s)).ToArray(),
                data.Item3.Select(s => Flatten(s)).ToArray(),
                data.Item4.Select(s => Flatten(s)).ToArray(),
            };

        internal static object[] Flatten(HTerm term) =>
            new object[]
            {
                term.Item1.ToArray(),
                term.Item2.ToArray()
            };

        internal static object[] Flatten(JordanWigner.JordanWignerInputState data) =>
            new object[]
            {
                new object[]
                {
                    data.Item1.Item1,
                    data.Item1.Item2
                },
                data.Item2.ToArray()
            };
    }

}
