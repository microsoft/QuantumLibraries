// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.IO;
using System.Text;
using System.Text.RegularExpressions;
using YamlDotNet.Core;
using YamlDotNet.Serialization;

namespace Microsoft.Quantum.Chemistry
{

    public partial class FermionHamiltonian
    {
        /// <summary>
        ///      Loads a Hamiltonian from integral data represented
        ///      in .yaml format.
        ///      Please see the <a href="http://ToBeDetermined/qchem-0.1.schema.json">
        ///      for further details about the
        ///      format parsed by this method.
        /// </summary>
        /// <param name="filename">The name of the file to be loaded.</param>
        /// <returns>
        ///      An instance of <see cref="FermionHamiltonian"/> representing the
        ///      data contained in <paramref name="filename"/>.
        /// </returns>
        public static IEnumerable<FermionHamiltonian> LoadFromYAML(string filename)
        {
            using (var reader = File.OpenText(filename))
            {
                var deserializer = new DeserializerBuilder().Build();
                var yamlData = deserializer.Deserialize<IntegralDataSchema>(reader);

                return LoadData.LoadIntegralData(yamlData);
            }
        }

    }


    internal struct Format
    {
        [YamlMember(Alias = "version", ApplyNamingConventions = false)]
        public string Version { get; set; }
    }

    internal struct Generator
    {
        [YamlMember(Alias = "source", ApplyNamingConventions = false)]
        public string Source { get; set; }

        [YamlMember(Alias = "version", ApplyNamingConventions = false)]
        public string Version { get; set; }
    }

    internal struct BasisSet
    {
        [YamlMember(Alias = "type", ApplyNamingConventions = false)]
        public string Type { get; set; }

        [YamlMember(Alias = "name", ApplyNamingConventions = false)]
        public string Name { get; set; }
    }

    public enum BibliographyKind
    {
        arXiv, DOI, URL
    }

    internal struct BibliographyItem : IYamlConvertible
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

    internal class HasUnits
    {
        [YamlMember(Alias = "units", ApplyNamingConventions = false)]
        // FIXME: make this an enum of allowed units.
        public string Units { get; set; }
    }

    internal class Geometry : HasUnits
    {
        [YamlMember(Alias = "coordinate_system", ApplyNamingConventions = false)]
        public string CoordinateSystem { get; set; }

        [YamlMember(Alias = "symmetry", ApplyNamingConventions = false)]
        public string Symmetry { get; set; }

        [YamlMember(Alias = "atoms", ApplyNamingConventions = false)]
        // FIXME: new struct for atoms
        public List<Dictionary<string, object>> Atoms { get; set; }
    }

    internal class SimpleQuantity : HasUnits
    {
        [YamlMember(Alias = "value", ApplyNamingConventions = false)]
        public double Value { get; set; }
    }

    internal class BoundedQuantity : HasUnits
    {
        [YamlMember(Alias = "value", ApplyNamingConventions = false)]
        public double? Value { get; set; }

        [YamlMember(Alias = "upper", ApplyNamingConventions = false)]
        public double Upper { get; set; }

        [YamlMember(Alias = "lower", ApplyNamingConventions = false)]
        public double Lower { get; set; }
    }

    internal class ArrayQuantity<TIndex, TValue> : HasUnits, IYamlConvertible
    {
        // TODO: make this an enum.
        public string Format { get; set; }
        public Dictionary<TIndex[], TValue> Values { get; set; }

        public void Read(IParser parser, Type expectedType, ObjectDeserializer nestedObjectDeserializer)
        {
            // We need to read units ourselves since we've overriden IYamlConvertible interface.
            var data = (Dictionary<string, object>)nestedObjectDeserializer(typeof(Dictionary<string, object>));
            Units = (string)data["units"];
            Format = (string)data["format"];
            Values = ((IEnumerable<object>)data["values"])
                .Select(entry => (IEnumerable<object>)entry)
                .ToDictionary(
                    entry => entry
                        .Take(entry.Count() - 1)
                        .Select((idx) => (TIndex)Convert.ChangeType(
                            idx,
                            typeof(TIndex)
                        ))
                        .ToArray(),
                    entry => (TValue)Convert.ChangeType(
                        entry.ElementAt(entry.Count() - 1),
                        typeof(TValue)
                    )
                );
        }

        public void Write(IEmitter emitter, ObjectSerializer nestedObjectSerializer)
        {
            nestedObjectSerializer(new Dictionary<string, object>
            {
                ["units"] = Units,
                ["format"] = Format,
                ["values"] = Values
                    .Select(kvp =>
                        kvp.Key.Select((idx) => (object)idx).Concat(new object[] { kvp.Value })
                    )
                    .ToList()
            });
        }
    }

    internal struct HamiltonianData
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

    internal struct IntegralSet
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

    internal struct SuggestedState
    {
        [YamlMember(Alias = "state", ApplyNamingConventions = false)]
        public SuggestedStateData SuggestedStateData { get; set; }
    }

    internal struct SuggestedStateData
    {
        [YamlMember(Alias = "label", ApplyNamingConventions = false)]
        public string Label { get; set; }

        [YamlMember(Alias = "energy", ApplyNamingConventions = false)]
        public SimpleQuantity Energy { get; set; }

        [YamlMember(Alias = "superposition", ApplyNamingConventions = false)]
        public List<List<object>> Superposition { get; set; }
    }

    internal struct IntegralDataSchema
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

    public partial class LoadData
    {
        internal static IEnumerable<FermionHamiltonian> LoadIntegralData(IntegralDataSchema schemaInstance, Double threshold = 1e-8)
        {
            return schemaInstance.IntegralSets.Select(
                (hamiltonianData, index) =>
                {
                    var hamiltonian = new FermionHamiltonian()
                    {

                        NOrbitals = hamiltonianData.NOrbitals,
                        EnergyOffset = 0 * hamiltonianData.EnergyOffset.Value + hamiltonianData.CoulombRepulsion.Value,
                        // TODO: optionally look at metadata to decide a better label.
                        Name = $"hamiltonian_{index}"
                    };

                    var fileHIJTerms = new Dictionary<Int64[], Double>(new Extensions.IntArrayIEqualityComparer());
                    var fileHIJKLTerms = new Dictionary<Int64[], Double>(new Extensions.IntArrayIEqualityComparer());

                    fileHIJTerms = hamiltonianData.Hamiltonian.OneElectronIntegrals.Values;
                    fileHIJKLTerms = hamiltonianData.Hamiltonian.TwoElectronIntegrals.Values;
                    hamiltonian.NOrbitals = hamiltonianData.NOrbitals;
                    hamiltonian.NElectrons = hamiltonianData.NElectrons;
                    HashSet<Int64[]> recordedIJTerms = new HashSet<Int64[]>(new Extensions.IntArrayIEqualityComparer());
                    foreach (var ijTerm in fileHIJTerms)
                    {
                        // If equivalent overlaps are specified, only record the first instance.
                        var orderedTerm = ijTerm.Key[0] > ijTerm.Key[1] ? ijTerm.Key.Reverse().ToArray() : ijTerm.Key;
                        if (recordedIJTerms.Contains(orderedTerm) == false)
                        {
                            recordedIJTerms.Add(orderedTerm);
                            if (Math.Abs(ijTerm.Value) >= threshold)
                            {
                                hamiltonian.AddFermionTerm(new OrbitalIntegral(ijTerm.Key.Select(o => o - 1), ijTerm.Value, OrbitalIntegral.Convention.Mulliken));
                            }
                        }
                    }
                    foreach (var ijklTerm in fileHIJKLTerms)
                    {
                        var tmp = new OrbitalIntegral(ijklTerm.Key.Select(o => o - 1), ijklTerm.Value, OrbitalIntegral.Convention.Mulliken);
                        if (Math.Abs(ijklTerm.Value) >= threshold)
                        {
                            hamiltonian.AddFermionTerm(new OrbitalIntegral(ijklTerm.Key.Select(o => o - 1), ijklTerm.Value, OrbitalIntegral.Convention.Mulliken));
                        }
                    }
                    hamiltonian.SortAndAccumulate();

                    if (!(hamiltonianData.SuggestedState == null))
                    {

                        foreach (var state in hamiltonianData.SuggestedState)
                        {
                            var label = state.SuggestedStateData.Label;
                            var energy = state.SuggestedStateData.Energy?.Value ?? 0.0;
                            var superpositionRaw = state.SuggestedStateData.Superposition;
                            var stringData = superpositionRaw.Select(o => o.Select(k => k.ToString()).ToList());
                            var superposition = stringData.Select(o => ParseInputState(o)).ToArray();
                            //Only have terms with non-zero amplitudes.
                            superposition = superposition.Where(o => Math.Abs(o.Item1.Item1) >= threshold).ToArray();
                            if (superposition.Count() > 0)
                            {
                                hamiltonian.InputStates.Add(
                                    new FermionHamiltonian.InputState { Label = label, Energy = energy, Superposition = superposition }
                                    );
                            }
                        }
                    }
                    return hamiltonian;
                }
            );
        }



        public static ((Double, Double), FermionTerm) ParseInputState(List<string> superpositionElement)
        {
            var amplitude = Double.Parse(superpositionElement.First(), System.Globalization.CultureInfo.InvariantCulture);
            var initialState = superpositionElement.Last();
            var ca = new List<Int64>();
            var so = new List<SpinOrbital>();

            for (int i = 1; i < superpositionElement.Count() - 1; i++)
            {
                FermionTerm singleTerm = ParsePolishNotation(superpositionElement[i]);
                ca.Add(singleTerm.CreationAnnihilationIndices.First());
                so.Add(singleTerm.SpinOrbitalIndices.First());
            }
            FermionTerm term = new FermionTerm(ca.ToArray(), so.ToArray(), amplitude);

            var canonicalOrder = term.ToCanonicalOrder();
            var noAnnihilationTerm = canonicalOrder.Where(o => !(o.CreationAnnihilationIndices.Contains(0)));
            var finalAmplitude = 0.0;
            // check if there are any valid terms.
            if (noAnnihilationTerm.Count() == 1)
            {
                term = noAnnihilationTerm.Single();
                finalAmplitude = term.coeff;
            }
            term.coeff = 1.0;
            return ((finalAmplitude, 0.0), term);
        }

        public static FermionTerm ParsePolishNotation(string input)
        {
            // Regex match examples: (1a)+ (2a)+ (3a)+ (4a)+ (5a)+ (6a)+ (1b)+ (2b)- (3b)+
            Regex regex = new Regex(@"(\((?<orbital>\d+)(?<spin>[ab])\)(?<operator>\+*))");
            Match match = regex.Match(input);
            if (match.Success)
            {
                var orbital = Int64.Parse(match.Groups["orbital"].ToString()) - 1;
                var spin = match.Groups["spin"].ToString() == "a" ? Spin.u : Spin.d;
                var conjugate = match.Groups["operator"].ToString() == "+" ? 1 : 0;
                return new FermionTerm(new long[] { conjugate }, new SpinOrbital[] { new SpinOrbital(orbital, spin) }, 1.0);
            }
            else
            {
                throw new System.ArgumentException($"{input} is not valid Polish notation");
            }
        }
    }
}