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

namespace Microsoft.Quantum.Chemistry
{
    /// <summary>
    /// Functionality to validate and parse Broombridge
    /// </summary>
    public static class Broombridge
    {
        /// <summary>
        /// Enumerable item for Broombridge version numbers.
        /// </summary>
        public enum Version
        {
            v0_1 = 0, v0_2 = 1
        }

        // Implement backwards compatibility by converting v0.1 to v0.2.
        // Only parse v0.2.

        /// <summary>
        /// Broombridge deserializers
        /// </summary>
        public static class Deserialize
        {
            /// <summary>
            /// Broombridge v0.1 deserializer
            /// </summary>
            /// <param name="filename">Broombridge filename to deserialize</param>
            /// <returns>Deserialized Broombridge v0.1 data.</returns>
            public static V0_1.Data v0_1(string filename)
            {
                using (var reader = File.OpenText(filename))
                {
                    var deserializer = new DeserializerBuilder().Build();
                    return deserializer.Deserialize<V0_1.Data>(reader);
                }
            }
            
            /// <summary>
            /// Broombridge deserializer v0.2
            /// </summary>
            /// <param name="filename">Broombridge filename to deserialize</param>
            /// <returns>Deserialized Broombridge v0.2 data.</returns>
            public static V0_2.Data v0_2(string filename)
            {
                using (var reader = File.OpenText(filename))
                {
                    var deserializer = new DeserializerBuilder().Build();
                    return deserializer.Deserialize<V0_2.Data>(reader);
                }
            }
        }

        /// <summary>
        /// Broombridge serializers
        /// </summary>
        public static class Serialize
        {
            /// <summary>
            /// Broombridge serializer
            /// </summary>
            /// <param name="filename">Broombridge filename to serialize</param>
            /// <returns>Serialized Broombridge</returns>
            public static void v0_2(V0_2.Data data, string filename)
            {
                var stringBuilder = new StringBuilder();
                var serializer = new Serializer();
                stringBuilder.AppendLine(serializer.Serialize(data));
                Console.WriteLine(stringBuilder);
                Console.WriteLine("");
            }
        }



        public static class Update
        {
            /// <summary>
            /// Converts v0.1 Broombridge to v0.2.
            /// </summary>
            /// <param name="input">Source Broombridge in v0.1 format.</param>
            /// <returns>Converted Broombridge in v0.2 format.</returns>
            public static V0_2.Data Data(V0_1.Data input)
            {
                var output = new V0_2.Data();

                output.Schema = input.Schema;
                output.Version = V0_2.Strings.VersionNumber;
                output.Generator = input.Generator;
                output.Bibliography = input.Bibliography;
                output.ProblemDescription = new List<V0_2.ProblemDescription>();

                foreach(var integralSet in input.IntegralSets)
                {
                    var problemDescription = new V0_2.ProblemDescription();
                    problemDescription.Metadata = integralSet.Metadata;
                    problemDescription.BasisSet = integralSet.BasisSet;
                    problemDescription.Geometry = integralSet.Geometry;
                    problemDescription.CoulombRepulsion = integralSet.CoulombRepulsion;
                    problemDescription.ScfEnergy = integralSet.ScfEnergy;
                    problemDescription.ScfEnergyOffset = integralSet.ScfEnergyOffset;
                    problemDescription.FciEnergy = integralSet.FciEnergy;
                    problemDescription.NOrbitals = integralSet.NOrbitals;
                    problemDescription.NElectrons = integralSet.NElectrons;
                    problemDescription.EnergyOffset = integralSet.EnergyOffset;
                    problemDescription.Hamiltonian = integralSet.Hamiltonian;

                    problemDescription.InitialStates = new List<V0_2.State>();
                    foreach(var sourceInitialState in integralSet.SuggestedState)
                    {
                        var initialState = new V0_2.State();
                        initialState.Label = sourceInitialState.SuggestedStateData.Label;
                        initialState.Energy = sourceInitialState.SuggestedStateData.Energy;
                        initialState.Method = V0_2.Strings.SparseMultiConfigurational;
                        initialState.Superposition = sourceInitialState.SuggestedStateData.Superposition;

                        problemDescription.InitialStates.Add(initialState);
                    }
                    
                    output.ProblemDescription.Add(problemDescription);
                }

                return output;
            }
        }

        public static class DataStructures
        {

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
                public DataStructures.ArrayQuantity<object, object> ParticleHoleRepresentation { get; set; }

                [YamlMember(Alias = "one_electron_integrals", ApplyNamingConventions = false)]
                public DataStructures.ArrayQuantity<long, double> OneElectronIntegrals { get; set; }

                [YamlMember(Alias = "two_electron_integrals", ApplyNamingConventions = false)]
                public DataStructures.ArrayQuantity<long, double> TwoElectronIntegrals { get; set; }

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

            /// <summary>
            /// Coupled-cluster operator
            /// </summary>
            public struct ClusterOperator
            {
                [YamlMember(Alias = "reference_state", ApplyNamingConventions = false)]
                public string Reference { get; set; }

                [YamlMember(Alias = "one_body_amplitudes", ApplyNamingConventions = false)]
                public List<List<string>> OneBodyAmplitudes { get; set; }

                [YamlMember(Alias = "two_body_amplitudes", ApplyNamingConventions = false)]
                public List<List<string>> TwoBodyAmplitudes { get; set; }
            }

        }

        /// <summary>
        /// Broombridge v0.2 format.
        /// 
        /// Changes from v0.1:
        /// - `format -> version` replaced with `version`.
        /// - `integral_sets` replaced with `problem_description`.
        /// - `metadata` is now optional.
        /// - Initial state suggestion `state` key removed. All values in this key are moved up one level.
        /// </summary>
        #region Broombridge v0.2 format
        public static class V0_2
        {
            public static class Strings {
                public static string SparseMultiConfigurational = "sparse_multi_configurational";
                public static string UnitaryCoupledCluster = "unitary_coupled_cluster";
                public static string VersionNumber = "0.2";
            }
            public static Version VersionNumber = Version.v0_2;

            // Root of Broombridge data structure
            public struct Data
            {

                [YamlMember(Alias = "$schema", ApplyNamingConventions = false)]
                public string Schema { get; set; }

                [YamlMember(Alias = "version", ApplyNamingConventions = false)]
                public string Version { get; set; }

                [YamlMember(Alias = "generator", ApplyNamingConventions = false)]
                public DataStructures.Generator Generator { get; set; }

                [YamlMember(Alias = "bibliography", ApplyNamingConventions = false)]
                public List<DataStructures.BibliographyItem> Bibliography { get; set; }

                [YamlMember(Alias = "problem_description", ApplyNamingConventions = false)]
                public List<ProblemDescription> ProblemDescription { get; set; }




            }

            public struct ProblemDescription
            {
                [YamlMember(Alias = "metadata", ApplyNamingConventions = false)]
                public Dictionary<string, object> Metadata { get; set; }

                [YamlMember(Alias = "basis_set", ApplyNamingConventions = false)]
                public DataStructures.BasisSet BasisSet { get; set; }

                [YamlMember(Alias = "geometry", ApplyNamingConventions = false)]
                public DataStructures.Geometry Geometry { get; set; }

                [YamlMember(Alias = "coulomb_repulsion", ApplyNamingConventions = false)]
                public DataStructures.SimpleQuantity CoulombRepulsion { get; set; }

                [YamlMember(Alias = "scf_energy", ApplyNamingConventions = false)]
                public DataStructures.SimpleQuantity ScfEnergy { get; set; }

                [YamlMember(Alias = "scf_energy_offset", ApplyNamingConventions = false)]
                public DataStructures.SimpleQuantity ScfEnergyOffset { get; set; }

                [YamlMember(Alias = "fci_energy", ApplyNamingConventions = false)]
                public DataStructures.BoundedQuantity FciEnergy { get; set; }

                [YamlMember(Alias = "n_orbitals", ApplyNamingConventions = false)]
                public int NOrbitals { get; set; }

                [YamlMember(Alias = "n_electrons", ApplyNamingConventions = false)]
                public int NElectrons { get; set; }

                [YamlMember(Alias = "energy_offset", ApplyNamingConventions = false)]
                public DataStructures.SimpleQuantity EnergyOffset { get; set; }

                [YamlMember(Alias = "hamiltonian", ApplyNamingConventions = false)]
                public DataStructures.HamiltonianData Hamiltonian { get; set; }

                // FIXME: actually specify what initial_state_suggestions looks like.
                //[YamlMember(Alias = "initial_state_suggestions", ApplyNamingConventions = false)]
                //public List<Dictionary<string, object>> InitialStateSuggestions { get; set; }

                [YamlMember(Alias = "initial_state_suggestions", ApplyNamingConventions = false)]
                public List<State> InitialStates { get; set; }
            }



            public struct State
            {
                [YamlMember(Alias = "label", ApplyNamingConventions = false)]
                public string Label { get; set; }

                [YamlMember(Alias = "method", ApplyNamingConventions = false)]
                public string Method { get; set; }

                [YamlMember(Alias = "energy", ApplyNamingConventions = false)]
                public DataStructures.SimpleQuantity Energy { get; set; }

                /// <summary>
                /// Sparse multi-configurational data
                /// </summary>
                [YamlMember(Alias = "superposition", ApplyNamingConventions = false)]
                public List<List<object>> Superposition { get; set; }

                /// <summary>
                /// Coupled-cluster operator
                /// </summary>
                [YamlMember(Alias = "cluster_operator", ApplyNamingConventions = false)]
                public DataStructures.ClusterOperator ClusterOperator { get; set; }
            }


        }
        #endregion

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
                public DataStructures.Generator Generator { get; set; }

                [YamlMember(Alias = "bibliography", ApplyNamingConventions = false)]
                public List<DataStructures.BibliographyItem> Bibliography { get; set; }


                [YamlMember(Alias = "integral_sets", ApplyNamingConventions = false)]
                public List<IntegralSet> IntegralSets { get; set; }

            }

            public struct Format
            {
                [YamlMember(Alias = "version", ApplyNamingConventions = false)]
                public string Version { get; set; }
            }


            public struct IntegralSet
            {
                [YamlMember(Alias = "metadata", ApplyNamingConventions = false)]
                public Dictionary<string, object> Metadata { get; set; }

                [YamlMember(Alias = "basis_set", ApplyNamingConventions = false)]
                public DataStructures.BasisSet BasisSet { get; set; }

                [YamlMember(Alias = "geometry", ApplyNamingConventions = false)]
                public DataStructures.Geometry Geometry { get; set; }

                [YamlMember(Alias = "coulomb_repulsion", ApplyNamingConventions = false)]
                public DataStructures.SimpleQuantity CoulombRepulsion { get; set; }

                [YamlMember(Alias = "scf_energy", ApplyNamingConventions = false)]
                public DataStructures.SimpleQuantity ScfEnergy { get; set; }

                [YamlMember(Alias = "scf_energy_offset", ApplyNamingConventions = false)]
                public DataStructures.SimpleQuantity ScfEnergyOffset { get; set; }

                [YamlMember(Alias = "fci_energy", ApplyNamingConventions = false)]
                public DataStructures.BoundedQuantity FciEnergy { get; set; }

                [YamlMember(Alias = "n_orbitals", ApplyNamingConventions = false)]
                public int NOrbitals { get; set; }

                [YamlMember(Alias = "n_electrons", ApplyNamingConventions = false)]
                public int NElectrons { get; set; }

                [YamlMember(Alias = "energy_offset", ApplyNamingConventions = false)]
                public DataStructures.SimpleQuantity EnergyOffset { get; set; }

                [YamlMember(Alias = "hamiltonian", ApplyNamingConventions = false)]
                public DataStructures.HamiltonianData Hamiltonian { get; set; }

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
                public DataStructures.SimpleQuantity Energy { get; set; }

                [YamlMember(Alias = "superposition", ApplyNamingConventions = false)]
                public List<List<object>> Superposition { get; set; }
            }


            #endregion


        }


        public struct ProblemDescription
        {


            /* public bool Energy_Provided;
             public double Energy_Min;
             public double Energy_Approx;
             public double Energy_Max;
             */
            public int NOrbitals, NElectrons;

            public double IdentityTerm;
            public List<OrbitalIntegral> OneBodyTerms;
            public List<OrbitalIntegral> TwoBodyTerms;
            /*[YamlMember(Alias = "coulomb_repulsion", ApplyNamingConventions = false)]
            public DataStructures.SimpleQuantity CoulombRepulsion { get; set; }

            [YamlMember(Alias = "scf_energy", ApplyNamingConventions = false)]
            public DataStructures.SimpleQuantity ScfEnergy { get; set; }

            [YamlMember(Alias = "scf_energy_offset", ApplyNamingConventions = false)]
            public DataStructures.SimpleQuantity ScfEnergyOffset { get; set; }

            [YamlMember(Alias = "fci_energy", ApplyNamingConventions = false)]
            public DataStructures.BoundedQuantity FciEnergy { get; set; }

            [YamlMember(Alias = "n_orbitals", ApplyNamingConventions = false)]
            public int NOrbitals { get; set; }

            [YamlMember(Alias = "n_electrons", ApplyNamingConventions = false)]
            public int NElectrons { get; set; }

            [YamlMember(Alias = "energy_offset", ApplyNamingConventions = false)]
            public DataStructures.SimpleQuantity EnergyOffset { get; set; }

            [YamlMember(Alias = "hamiltonian", ApplyNamingConventions = false)]
            public DataStructures.HamiltonianData Hamiltonian { get; set; }
    */
            // FIXME: actually specify what initial_state_suggestions looks like.
            //[YamlMember(Alias = "initial_state_suggestions", ApplyNamingConventions = false)]
            //public List<Dictionary<string, object>> InitialStateSuggestions { get; set; }

            //[YamlMember(Alias = "initial_state_suggestions", ApplyNamingConventions = false)]
            public ProblemDescription(Broombridge.V0_2.ProblemDescription broombridgeProblem)
            {
                IdentityTerm = broombridgeProblem.CoulombRepulsion.Value + broombridgeProblem.EnergyOffset.Value;

                NOrbitals = broombridgeProblem.NOrbitals;
                NElectrons = broombridgeProblem.NElectrons;

                var test = broombridgeProblem.Hamiltonian.OneElectronIntegrals.Values
    
        }
        }
    }


    public static partial class Extensions
    {
        /// <summary>
        ///      Given a value of an enumeration type, and an action for each
        ///      possible value of that enumeration type, performs the action
        ///      corresponding to the given value.
        /// </summary>
        public static void Map<E>
            (this E @enum, params (E, Action)[] actions)
            where E : struct, IConvertible
        {
            foreach (var (value, action) in actions)
            {
                if (@enum.Equals(value))
                {
                    action();
                    return;
                }
            }
        }

        /// <summary>
        ///      Given a value of an enumeration type, and an function for each
        ///      possible value of that enumeration type, returns the value
        //       returned by the function corresponding to the given value.
        /// </summary>
        public static T Map<E, T>
            (this E @enum, params (E, Func<T>)[] actions)
            where E : struct, IConvertible
        {
            foreach (var (value, action) in actions)
            {
                if (@enum.Equals(value))
                {
                    return action();
                }
            }

            throw new ArgumentException($"Expected {@enum} to be a member of {@enum.GetType()}.");
        }
    }
}