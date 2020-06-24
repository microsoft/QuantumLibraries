// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

#nullable enable

using Microsoft.Quantum.Simulation.Core;

using System.Text.RegularExpressions;
using System.Linq;
using System.Collections.Generic;

using Microsoft.Quantum.Chemistry.OrbitalIntegrals;
using System.IO;
using System;

namespace Microsoft.Quantum.Chemistry
{


    public static class LiQuiDSerializer
    {
        public static IEnumerable<ElectronicStructureProblem> Deserialize(TextReader reader)
        {
            var allText = reader.ReadToEnd();
            string[] delimiters = { "tst" };
            var lines = allText.Split(delimiters, System.StringSplitOptions.RemoveEmptyEntries);
            var hamiltonians = lines.Select(DeserializeSingleProblem);
            return hamiltonians;
        }

        internal static ElectronicStructureProblem DeserializeSingleProblem(string line)
        {
            var problem = new ElectronicStructureProblem()
            {
                Metadata = new Dictionary<string, object>()
            };

            var regexMiscellaneous = new Regex(@"((info=(?<info>[^\s]*)))");
            var regexnuc = new Regex(@"nuc=(?<nuc>-?\s*\d*.\d*)");
            var regexPQ = new Regex(@"(^|\s+)(?<p>\d+),(?<q>\d+)\D*=\s*(?<coeff>-?\s*\d*.\d*)e?(?<exponent>-?\d*)");
            var regexPQRS = new Regex(@"(^|\s+)(?<p>\d+)\D+(?<q>\d+)\D+(?<r>\d+)\D+(?<s>\d+)\D*=\s*(?<coeff>-?\s*\d*.\d*)e?(?<exponent>-?\d*)");

            double coulombRepulsion = 0.0;
            var fileHIJTerms = new Dictionary<int[], double>(new Extensions.ArrayEqualityComparer<int>());
            var fileHIJKLTerms = new Dictionary<int[], double>(new Extensions.ArrayEqualityComparer<int>());
            var hamiltonian = new OrbitalIntegralHamiltonian();

            var nOrbitals = 0L;

            Match stringMisc = regexMiscellaneous.Match(line);
            if (stringMisc.Success)
            {
                problem.Metadata["misc_info"] = stringMisc.Groups["info"].ToString();
            }


            Match stringnuc = regexnuc.Match(line);
            if (stringnuc.Success)
            {
                coulombRepulsion = double.Parse(stringnuc.Groups["nuc"].ToString()).ToDoubleCoeff();
                hamiltonian.Add(TermType.OrbitalIntegral.Identity, new OrbitalIntegral(),  coulombRepulsion);
            }
            foreach (Match stringPQ in regexPQ.Matches(line))
            {
                if (stringPQ.Success)
                {
                    var p = int.Parse(stringPQ.Groups["p"].ToString());
                    var q = int.Parse(stringPQ.Groups["q"].ToString());
                    var coeff = double.Parse(stringPQ.Groups["coeff"].ToString());
                    var exponentString = stringPQ.Groups["exponent"].ToString();
                    var exponent = 0.0;
                    if (exponentString != "")
                    {
                        exponent = double.Parse(stringPQ.Groups["exponent"].ToString());
                    }
                    nOrbitals = new long[] { nOrbitals, p + 1, q + 1 }.Max();

                    var orbitalIntegral = new OrbitalIntegral(new int[] { p, q }, coeff * (10.0).Pow(exponent));
                    var orbitalIntegralCanonical = orbitalIntegral.ToCanonicalForm();
                    

                    if (fileHIJTerms.ContainsKey(orbitalIntegralCanonical.OrbitalIndices))
                    {
                        // Check consistency
                        if (fileHIJTerms[orbitalIntegralCanonical.OrbitalIndices] != orbitalIntegral.Coefficient)
                        {
                            // Consistency check failed.
                            throw new System.NotSupportedException(
                                $"fileHPQTerm Consistency check fail. " +
                                $"Orbital integral {orbitalIntegral} coefficient {orbitalIntegral.Coefficient}" +
                                $"does not match recorded {orbitalIntegralCanonical} coefficient {fileHIJTerms[orbitalIntegral.OrbitalIndices]}.");
                        }
                        else
                        {
                            // Consistency check passed.
                        }
                        
                    }
                    else
                    {
                        fileHIJTerms.Add(orbitalIntegralCanonical.OrbitalIndices, orbitalIntegralCanonical.Coefficient);
                    }
                }
            }
            foreach (Match stringPQRS in regexPQRS.Matches(line))
            {
                if (stringPQRS.Success)
                {
                    var p = int.Parse(stringPQRS.Groups["p"].ToString());
                    var q = int.Parse(stringPQRS.Groups["q"].ToString());
                    var r = int.Parse(stringPQRS.Groups["r"].ToString());
                    var s = int.Parse(stringPQRS.Groups["s"].ToString());
                    var coeff = double.Parse(stringPQRS.Groups["coeff"].ToString());
                    var exponentString = stringPQRS.Groups["exponent"].ToString();
                    var exponent = 0.0;
                    if (exponentString != "")
                    {
                        exponent = double.Parse(stringPQRS.Groups["exponent"].ToString());
                    }
                    nOrbitals = new long[] { nOrbitals, p + 1, q + 1, r + 1, s + 1 }.Max();

                    var orbitalIntegral = new OrbitalIntegral(new int[] { p, q, r, s }, coeff * (10.0).Pow(exponent));
                    var orbitalIntegralCanonical = orbitalIntegral.ToCanonicalForm();
                    

                    if (fileHIJKLTerms.ContainsKey(orbitalIntegralCanonical.OrbitalIndices))
                    {
                        // Check consistency
                        if (fileHIJKLTerms[orbitalIntegralCanonical.OrbitalIndices] != orbitalIntegral.Coefficient)
                        {
                            // Consistency check failed.
                            throw new System.NotSupportedException(
                                $"fileHPQRSTerm Consistency check fail. " +
                                $"Orbital integral {orbitalIntegral.OrbitalIndices} coefficient {orbitalIntegral.Coefficient}" +
                                $"does not match recorded {orbitalIntegralCanonical.OrbitalIndices} coefficient {fileHIJKLTerms[orbitalIntegral.OrbitalIndices]}.");
                        }
                        else
                        {
                            // Consistency check passed.
                        }
                    }
                    else
                    {
                        fileHIJKLTerms.Add(orbitalIntegralCanonical.OrbitalIndices, orbitalIntegralCanonical.Coefficient);
                    }
                }
            }

            hamiltonian.Add(fileHIJTerms.Select(o => new OrbitalIntegral(o.Key, o.Value)).ToList());

            hamiltonian.Add(fileHIJKLTerms.Select(o => new OrbitalIntegral(o.Key, o.Value)).ToList());

            problem.OrbitalIntegralHamiltonian = hamiltonian;
            problem.NOrbitals = System.Convert.ToInt32 ( nOrbitals);
            problem.CoulombRepulsion = coulombRepulsion.WithUnits("hartree");

            return problem;
        }

        public static void Serialize(TextWriter writer, IEnumerable<ElectronicStructureProblem> problems)
        {
            throw new NotImplementedException("Serialization to LiQuiD is not yet implemented.");
        }
    }



    /// <summary>
    /// Methods for loading Hamiltonian data from standard formats
    /// into a <see cref="FermionHamiltonian"/>.
    /// </summary>
    [Obsolete(
        "Please use LiQuiDSerializer instead.",
        false
    )]

    public class LiQuiD
    {
        public struct ProblemDescription
        {
            public int NOrbitals { get; set; }
            public int NElectrons { get; set; }
            public double CoulombRepulsion { get; set; }
            public string MiscellaneousInformation { get; set; }
            public OrbitalIntegralHamiltonian OrbitalIntegralHamiltonian { get; set; }
        }

        /// <summary>
        ///      Loads a Hamiltonian from integral data represented
        ///      in LIQ𝑈𝑖|⟩ format.
        ///      Please see the <a href="https://stationq.github.io/Liquid/docs/LIQUiD.pdf">
        ///      LIQ𝑈𝑖|⟩ documentation</a> for further details about the
        ///      format parsed by this method.
        /// </summary>
        /// <param name="filename">The name of the file to be loaded.</param>
        /// <returns>
        ///      List of electronic structure problem deserialized from the file.
        /// </returns>
        public static IEnumerable<ProblemDescription> Deserialize(string filename)
        {
            var name = filename;
            using var reader = File.OpenText(filename);
            return LiQuiDSerializer
                .Deserialize(reader)
                .Select(problem => new ProblemDescription
                {
                    CoulombRepulsion = problem.CoulombRepulsion.Value,
                    MiscellaneousInformation = problem.Metadata.GetValueOrDefault("misc_info", "").ToString(),
                    NElectrons = problem.NElectrons,
                    NOrbitals = problem.NOrbitals,
                    OrbitalIntegralHamiltonian = problem.OrbitalIntegralHamiltonian
                });
        }

    }
}
