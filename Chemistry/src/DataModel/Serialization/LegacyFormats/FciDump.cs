// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using Microsoft.Quantum.Simulation.Core;

using System.Text.RegularExpressions;
using System.Linq;
using System.Collections.Generic;
using System.IO;

using Microsoft.Quantum.Chemistry.OrbitalIntegrals;
using static Microsoft.Quantum.Chemistry.OrbitalIntegrals.IndexConventionConversions;
using System;
using System.Diagnostics;

namespace Microsoft.Quantum.Chemistry
{
    /// <summary>
    ///      Serialization and deserialization support for FCIDUMP
    ///      formatted problem descriptions.
    /// </summary>
    public static class FciDumpSerializer
    {

        /// <summary>
        ///     Deserializes an FCIDUMP-formatted problem description.
        /// </summary>
        /// <param name="reader">A stream for reading FCIDUMP data.</param>
        /// <returns>
        ///      An electronic structure problem deserialized from the file.
        /// </returns>
        public static IEnumerable<ElectronicStructureProblem> Deserialize(TextReader reader)
        {
            // FCIDUMP files begin with a FORTRAN-formatted namelist, delimited
            // by &FCI and &END. We start by extracting that namelist.
            var allText = reader.ReadToEnd();
            var lines = Regex.Split(allText, "\r\n|\r|\n");
            if (lines == null)
            {
                throw new IOException("Expected a non-empty FCIDUMP file.");
            }
            var header = System.String.Join("\n", lines.TakeWhile(line => line.Trim() != "&END")).Trim();
            var body = lines!.SkipWhile(line => line.Trim() != "&END").Skip(1).ToList();
            
            // Make sure that the header starts with &FCI, as expected.
            if (!header.StartsWith("&FCI"))
            {
                throw new IOException("FCIDUMP file did not start with \"&FCI\" as expected.");
            }
            
            // Split out the &FCI and &END lines, turn the rest into a dictionary of namelist items.
            var namelist = Regex.Matches(
                header
                .Replace("&FCI", "")
                .Replace("&END", ""),
                pattern: "\\s*(?<identifier>\\w+)\\s*=\\s*(?<value>[^=]+),\\s*"
            )
            .ToDictionary(
                match => match.Groups["identifier"].Value,
                match => match.Groups["value"].Value
            );

            var hamiltonian = new OrbitalIntegralHamiltonian();
            var arrayData = body
                .Select(line => line.Trim())
                .Where(line => line.Length > 0)
                .Select(
                    line => line.Split(" ", StringSplitOptions.RemoveEmptyEntries)
                )
                .Select(
                    row => (
                        Double.Parse(row[0]), 
                        row[1..].Select(Int32.Parse).Where(idx => idx != 0).ToZeroBasedIndices()
                    )
                );
            var (energyOffset, _) = arrayData.Where(item => item.Item2.Length == 0).Single();
            hamiltonian.Add(arrayData
                .Where(row => row.Item2.Length > 0)
                .SelectMaybe(
                    row => row.Item2.Length % 2 == 0
                           ? new OrbitalIntegral(
                                 row.Item2, row.Item1, OrbitalIntegral.Convention.Mulliken
                             ).ToCanonicalForm()
                           : null
                )
                .Distinct()
            );
            
            return new List<ElectronicStructureProblem>
            {
                new ElectronicStructureProblem
                {
                    EnergyOffset = energyOffset.WithUnits("hartree"),
                    Metadata = new Dictionary<string, object>
                    {
                        ["Comment"] = "Imported from FCIDUMP"
                    },
                    NElectrons = Int32.Parse(namelist["NELEC"]),
                    NOrbitals = Int32.Parse(namelist["NORB"]),
                    OrbitalIntegralHamiltonian = hamiltonian
                }
            };
        }

        public static void Serialize(TextWriter writer, IEnumerable<ElectronicStructureProblem> problems)
        {
            var problem = problems.Single();

            // Start by writing the header.
            writer.WriteLine($"&FCI NORB={problem.NOrbitals},NELEC={problem.NElectrons},");
            // Assume global phase symmetry for now.
            writer.WriteLine($" ORBSYM={String.Join("", Enumerable.Range(0, problem.NOrbitals).Select(idx => "1,"))}");
            writer.WriteLine($" ISYM=1,");
            writer.WriteLine("&END");

            (string, double) FormatTerm(KeyValuePair<OrbitalIntegral, DoubleCoeff> term) =>
                (
                    String.Join(" ",
                        ConvertIndices(
                            term.Key.OrbitalIndices,
                            OrbitalIntegral.Convention.Dirac,
                            OrbitalIntegral.Convention.Mulliken
                        )
                        .ToOneBasedIndices()
                        .Select(i => i.ToString())
                    ),
                    term.Value
                );

            foreach (var (indices, value) in problem
                                 .OrbitalIntegralHamiltonian
                                 .Terms[TermType.OrbitalIntegral.TwoBody]
                                 .Select(FormatTerm))
            {
                writer.WriteLine($"{value} {indices}");
            }
            // Next write out all one-body terms, using trailing zeros to indicate one-body.
            foreach (var (indices, value) in problem
                                 .OrbitalIntegralHamiltonian
                                 .Terms[TermType.OrbitalIntegral.OneBody]
                                 .Select(FormatTerm))
            {
                writer.WriteLine($"{value} {indices} 0 0");
            }
            // Finish by writing out the identity term.
            var identityTerm = problem
                .OrbitalIntegralHamiltonian
                .Terms.GetValueOrDefault(TermType.OrbitalIntegral.Identity, null)
                ?.Select(term => term.Value.Value)
                ?.SingleOrDefault() ?? 0.0
                + problem.EnergyOffset.Value;
            writer.WriteLine($"{identityTerm} 0 0 0 0");
        }
    }
}
