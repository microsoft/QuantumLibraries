// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

#nullable enable

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

using static Microsoft.Quantum.Chemistry.OrbitalIntegrals.IndexConventionConversions;

namespace Microsoft.Quantum.Chemistry.Broombridge
{

    internal static class BroombridgeExtensionsV0_1
    {

        internal static Broombridge.V0_1.SimpleQuantity ToBroombridgeV0_1(this Quantity<double> quantity) =>
            new Broombridge.V0_1.SimpleQuantity
            {
                Units = quantity.Units,
                Value = quantity.Value
            };

        internal static Broombridge.V0_1.BoundedQuantity ToBroombridgeV0_1(this BoundedQuantity<double> quantity) =>
            new Broombridge.V0_1.BoundedQuantity
            {
                Units = quantity.Units,
                Value = quantity.Value,
                Lower = quantity.Lower,
                Upper = quantity.Upper
            };

        internal static Quantity<double> FromBroombridgeV0_1(this Broombridge.V0_1.SimpleQuantity quantity) =>
            new Quantity<double>
            {
                Units = quantity.Units,
                Value = quantity.Value
            };

        internal static BoundedQuantity<double> FromBroombridgeV0_1(this Broombridge.V0_1.BoundedQuantity quantity) =>
            new BoundedQuantity<double>
            {
                Units = quantity.Units,
                Value = quantity.Value,
                Lower = quantity.Lower,
                Upper = quantity.Upper
            };

        internal static BasisSet FromBroombridgeV0_1(this Broombridge.V0_1.BasisSet basisSet) =>
            new BasisSet
            {
                Name = basisSet.Name,
                Type = basisSet.Type
            };

        

        internal static Geometry FromBroombridgeV0_1(this Broombridge.V0_1.Geometry geometry) =>
            new Geometry
            {
                Atoms = geometry.Atoms,
                CoordinateSystem = geometry.CoordinateSystem,
                Symmetry = geometry.Symmetry,
                Units = geometry.Units
            };

        internal static V0_1.Geometry ToBroombridgeV0_1(this Geometry geometry) =>
            new V0_1.Geometry
            {
                Atoms = geometry.Atoms,
                CoordinateSystem = geometry.CoordinateSystem,
                Symmetry = geometry.Symmetry,
                Units = geometry.Units
            };

        internal static V0_1.ArrayQuantity<long, double> ToBroombridgeV0_1(
                this Dictionary<OrbitalIntegrals.OrbitalIntegral, DoubleCoeff> terms
        ) =>
            new V0_1.ArrayQuantity<long, double>()
            {
                Format = "sparse",
                Units = "hartree",
                Values = terms.Select(term => (
                             ConvertIndices(
                                 term
                                .Key
                                .ToCanonicalForm()
                                .OrbitalIndices,
                                OrbitalIntegral.Convention.Dirac,
                                OrbitalIntegral.Convention.Mulliken
                             )
                             .ToOneBasedIndices()
                             .Select(idx => (long)idx)
                             .ToArray(),
                             term.Value.Value
                         )).ToList()
            };

        internal static V0_1.HamiltonianData ToBroombridgeV0_1(this OrbitalIntegralHamiltonian hamiltonian) =>
            new V0_1.HamiltonianData
            {
                OneElectronIntegrals = hamiltonian
                    .Terms[TermType.OrbitalIntegral.OneBody]
                    .ToBroombridgeV0_1(),
                TwoElectronIntegrals = hamiltonian
                    .Terms[TermType.OrbitalIntegral.TwoBody]
                    .ToBroombridgeV0_1()
            };

        internal static V0_1.IntegralSet ToBroombridgeV0_1(this ElectronicStructureProblem problem) =>
            throw new NotImplementedException("Not yet implemented.");

    }

    /// <summary>
    /// Broombridge v0.1 format
    /// </summary>
    #region Broombridge v0.1 format
    internal static class V0_1
    {
        public static string SchemaUrl = "https://raw.githubusercontent.com/microsoft/Quantum/master/Chemistry/Schema/broombridge-0.1.schema.json";

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
            public BasisSet? BasisSet { get; set; }

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

                Kind = kind.ToLowerInvariant() switch
                {
                    "arxiv" => BibliographyKind.arXiv,
                    "doi" => BibliographyKind.DOI,
                    _ => BibliographyKind.URL
                };

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

        public class BasisSet
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
            // TODO: Placeholder object for ParticleHoleRepresentation, which we do not
            // yet support.
            public ArrayQuantity<object, object>? ParticleHoleRepresentation { get; set; }

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
                    var a = entries.Take(entries.Count() - 1).Select(e => (TIndex)System.Convert.ChangeType(e, typeof(TIndex))).ToArray();
                    var q = (TValue)System.Convert.ChangeType(entries.Last(), typeof(TValue));
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
