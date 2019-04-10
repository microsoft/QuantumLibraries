// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Text;
using System.Linq;
using System.Text.RegularExpressions;
using YamlDotNet.Core;
using YamlDotNet.Serialization;

using Microsoft.Quantum.Chemistry.OrbitalIntegrals;
using Microsoft.Quantum.Chemistry.LadderOperators;

namespace Microsoft.Quantum.Chemistry.Broombridge
{
    /// <summary>
    /// Broombridge v0.1 format
    /// </summary>
    #region Broombridge v0.1 format
    public static class V0_1
    {
        public struct Data
        {

            [YamlMember(Alias = "$schema", ApplyNamingConventions = false)]
            public string Schema { get; set; }

            [YamlMember(Alias = "format", ApplyNamingConventions = false)]
            public Format Format { get; set; }

            [YamlMember(Alias = "generator", ApplyNamingConventions = false)]
            public Generator Generator { get; set; }

            [YamlMember(Alias = "bibliography", ApplyNamingConventions = false)]
            public List<BibliographyItem> Bibliography { get; set; }


            [YamlMember(Alias = "integral_sets", ApplyNamingConventions = false)]
            public List<IntegralSet> IntegralSets { get; set; }

        }



        public struct IntegralSet
        {
            [YamlMember(Alias = "metadata", ApplyNamingConventions = false)]
            public Dictionary<string, object> Metadata { get; set; }

            [YamlMember(Alias = "basis_set", ApplyNamingConventions = false)]
            public BasisSet BasisSet { get; set; }

            [YamlMember(Alias = "geometry", ApplyNamingConventions = false)]
            public Geometry Geometry { get; set; }

            [YamlMember(Alias = "coulomb_repulsion", ApplyNamingConventions = false)]
            public SimpleQuantity CoulombRepulsion { get; set; }

            [YamlMember(Alias = "scf_energy", ApplyNamingConventions = false)]
            public SimpleQuantity ScfEnergy { get; set; }

            [YamlMember(Alias = "scf_energy_offset", ApplyNamingConventions = false)]
            public SimpleQuantity ScfEnergyOffset { get; set; }

            [YamlMember(Alias = "fci_energy", ApplyNamingConventions = false)]
            public BoundedQuantity FciEnergy { get; set; }

            [YamlMember(Alias = "n_orbitals", ApplyNamingConventions = false)]
            public int NOrbitals { get; set; }

            [YamlMember(Alias = "n_electrons", ApplyNamingConventions = false)]
            public int NElectrons { get; set; }

            [YamlMember(Alias = "energy_offset", ApplyNamingConventions = false)]
            public SimpleQuantity EnergyOffset { get; set; }

            [YamlMember(Alias = "hamiltonian", ApplyNamingConventions = false)]
            public HamiltonianData Hamiltonian { get; set; }

            // FIXME: actually specify what initial_state_suggestions looks like.
            //[YamlMember(Alias = "initial_state_suggestions", ApplyNamingConventions = false)]
            //public List<Dictionary<string, object>> InitialStateSuggestions { get; set; }

            [YamlMember(Alias = "initial_state_suggestions", ApplyNamingConventions = false)]
            public List<SuggestedState> SuggestedState { get; set; }
        }

        public struct SuggestedState
        {
            [YamlMember(Alias = "state", ApplyNamingConventions = false)]
            public SuggestedStateData SuggestedStateData { get; set; }
        }

        public struct SuggestedStateData
        {
            [YamlMember(Alias = "label", ApplyNamingConventions = false)]
            public string Label { get; set; }

            [YamlMember(Alias = "energy", ApplyNamingConventions = false)]
            public SimpleQuantity Energy { get; set; }

            [YamlMember(Alias = "superposition", ApplyNamingConventions = false)]
            public List<List<string>> Superposition { get; set; }
        }


        public struct Format
        {
            [YamlMember(Alias = "version", ApplyNamingConventions = false)]
            public string Version { get; set; }
        }

        public struct Generator
        {
            [YamlMember(Alias = "source", ApplyNamingConventions = false)]
            public string Source { get; set; }

            [YamlMember(Alias = "version", ApplyNamingConventions = false)]
            public string Version { get; set; }
        }


        public enum BibliographyKind
        {
            arXiv, DOI, URL
        }
        public struct BibliographyItem : IYamlConvertible
        {
            public string Value { get; set; }
            public BibliographyKind Kind { get; set; }

            public void Read(IParser parser, Type expectedType, ObjectDeserializer nestedObjectDeserializer)
            {
                var dict = (Dictionary<string, string>)nestedObjectDeserializer(typeof(Dictionary<string, string>));
                // There should only be one item.
                Debug.Assert(dict.Count == 1);
                var kind = dict.Keys.ElementAt(0);

                Value = dict[kind];
                if (kind.ToLowerInvariant() == "arxiv")
                {
                    Kind = BibliographyKind.arXiv;
                }
                else if (kind.ToLowerInvariant() == "doi")
                {
                    Kind = BibliographyKind.DOI;
                }
                else
                {
                    Kind = BibliographyKind.URL;
                }
            }

            public void Write(IEmitter emitter, ObjectSerializer nestedObjectSerializer)
            {
                nestedObjectSerializer(new Dictionary<string, string>
                {
                    [Kind.Map(
                        (BibliographyKind.arXiv, () => "arXiv"),
                        (BibliographyKind.DOI, () => "doi"),
                        (BibliographyKind.URL, () => "url")
                    )] = Value
                });
            }
        }

        public struct BasisSet
        {
            [YamlMember(Alias = "type", ApplyNamingConventions = false)]
            public string Type { get; set; }

            [YamlMember(Alias = "name", ApplyNamingConventions = false)]
            public string Name { get; set; }
        }

        public class HasUnits
        {
            [YamlMember(Alias = "units", ApplyNamingConventions = false)]
            // FIXME: make this an enum of allowed units.
            public string Units { get; set; }
        }

        public class Geometry : HasUnits
        {
            [YamlMember(Alias = "coordinate_system", ApplyNamingConventions = false)]
            public string CoordinateSystem { get; set; }

            [YamlMember(Alias = "symmetry", ApplyNamingConventions = false)]
            public string Symmetry { get; set; }

            [YamlMember(Alias = "atoms", ApplyNamingConventions = false)]
            // FIXME: new struct for atoms
            public List<Dictionary<string, object>> Atoms { get; set; }
        }

        public struct HamiltonianData
        {
            [YamlMember(Alias = "particle_hole_representation", ApplyNamingConventions = false)]
            // TODO: make this not object
            // FIXME: currently strips off the last element as the "value", but the
            //        present schema requires us to pull off the last two as a (double, string).
            public ArrayQuantity<object, object> ParticleHoleRepresentation { get; set; }

            [YamlMember(Alias = "one_electron_integrals", ApplyNamingConventions = false)]
            public ArrayQuantity<long, double> OneElectronIntegrals { get; set; }

            [YamlMember(Alias = "two_electron_integrals", ApplyNamingConventions = false)]
            public ArrayQuantity<long, double> TwoElectronIntegrals { get; set; }

        }

        public class SimpleQuantity : HasUnits
        {
            [YamlMember(Alias = "value", ApplyNamingConventions = false)]
            public double Value { get; set; }
        }

        public class BoundedQuantity : HasUnits
        {
            [YamlMember(Alias = "value", ApplyNamingConventions = false)]
            public double? Value { get; set; }

            [YamlMember(Alias = "upper", ApplyNamingConventions = false)]
            public double Upper { get; set; }

            [YamlMember(Alias = "lower", ApplyNamingConventions = false)]
            public double Lower { get; set; }
        }

        public class ArrayQuantity<TIndex, TValue> : HasUnits, IYamlConvertible
        {
            // TODO: make this an enum.
            public string Format { get; set; }
            public List<(TIndex[], TValue)> Values { get; set; }

            public void Read(IParser parser, Type expectedType, ObjectDeserializer nestedObjectDeserializer)
            {
                // We need to read units ourselves since we've overriden IYamlConvertible interface.
                var data = (Dictionary<string, object>)nestedObjectDeserializer(typeof(Dictionary<string, object>));
                Units = (string)data["units"];
                Format = (string)data["format"];
                Values = new List<(TIndex[], TValue)>();
                foreach (var value in ((IEnumerable<object>)data["values"]))
                {
                    var entries = (IEnumerable<object>)value;
                    var a = entries.Take(entries.Count() - 1).Select(e => (TIndex)Convert.ChangeType(e, typeof(TIndex))).ToArray();
                    var q = (TValue)Convert.ChangeType(entries.Last(), typeof(TValue));
                    Values.Add((a, q));
                }
            }

            public void Write(IEmitter emitter, ObjectSerializer nestedObjectSerializer)
            {
                nestedObjectSerializer(new Dictionary<string, object>
                {
                    ["units"] = Units,
                    ["format"] = Format,
                    ["values"] = Values
                        .Select(entry =>
                            entry.Item1.Select((idx) => (object)idx).Concat(new object[] { entry.Item2 })
                        )
                        .ToList()
                });
            }
        }

        /// <summary>
        /// Builds Hamiltonian from Broombridge orbital integral data.
        /// </summary>
        internal static OrbitalIntegralHamiltonian ToOrbitalIntegralHamiltonian(
            HamiltonianData hamiltonianData)
        {
            var hamiltonian = new OrbitalIntegralHamiltonian();

            // This will convert from Broombridge 1-indexing to 0-indexing.
            hamiltonian.Add
                (hamiltonianData.OneElectronIntegrals.Values
                .Select(o => new OrbitalIntegral(o.Item1
                .Select(k => (int)(k - 1)), o.Item2, OrbitalIntegral.Convention.Mulliken)
                .ToCanonicalForm())
                .Distinct());

            // This will convert from Broombridge 1-indexing to 0-indexing.
            // This will convert to Dirac-indexing.
            hamiltonian.Add
                (hamiltonianData.TwoElectronIntegrals.Values
                .Select(o => new OrbitalIntegral(o.Item1
                .Select(k => (int)(k - 1)), o.Item2, OrbitalIntegral.Convention.Mulliken)
                .ToCanonicalForm())
                .Distinct());

            return hamiltonian;
        }

        // orbital, spin, raising/lowering.
        internal static LadderOperator<SpinOrbital> ParsePolishNotation(string input)
        {
            // Regex match examples: (1a)+ (2a)+ (3a)+ (4a)+ (5a)+ (6a)+ (1b)+ (2b)- (3b)+
            Regex regex = new Regex(@"(\((?<orbital>\d+)(?<spin>[ab])\)(?<operator>\+*))");
            Match match = regex.Match(input);
            if (match.Success)
            {
                // Convert from Broombridge 1-indexing to 0-indexing.
                var orbital = int.Parse(match.Groups["orbital"].ToString()) - 1;
                var spin = match.Groups["spin"].ToString() == "a" ? Spin.u : Spin.d;
                var conjugate = match.Groups["operator"].ToString() == "+" ? RaisingLowering.u : RaisingLowering.d;
                return new LadderOperator<SpinOrbital>(conjugate, new SpinOrbital(orbital, spin));
            }
            else
            {
                throw new System.ArgumentException($"{input} is not valid Polish notation");
            }
        }
    }

    #endregion
}