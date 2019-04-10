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
using System.Numerics;

using Microsoft.Quantum.Chemistry.OrbitalIntegrals;
using Microsoft.Quantum.Chemistry.Fermion;
using Microsoft.Quantum.Chemistry.Generic;

namespace Microsoft.Quantum.Chemistry.Broombridge
{
    // We enumerate these quantities to know which quantities
    // exactly match the previous version.
    using Format = V0_1.Format;
    using Generator = V0_1.Generator;
    using BibliographyItem = V0_1.BibliographyItem;
    using Geometry = V0_1.Geometry;
    using HasUnits = V0_1.HasUnits;
    using BasisSet = V0_1.BasisSet;
    using SimpleQuantity = V0_1.SimpleQuantity;
    using HamiltonianData = V0_1.HamiltonianData;
    using BoundedQuantity = V0_1.BoundedQuantity;


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
        // Lower-case invariant dictionary key comparer.
        private static Dictionary<string, StateType> StateTypeDictionary = new Dictionary<string, StateType>(StringComparer.OrdinalIgnoreCase)
        {
            { "single_configurational", StateType.SingleConfigurational },
            { "sparse_multi_configurational", StateType.SparseMultiConfigurational },
            { "unitary_coupled_cluster", StateType.SparseMultiConfigurational }
        };
            
        internal static class UpdaterStrings
        {
            public const string SingleConfigurational = "single_configurational";
            public const string SparseMultiConfigurational = "sparse_multi_configurational";
            public const string UnitaryCoupledCluster = "unitary_coupled_cluster";
            public const string VersionNumber = "0.2";
        }

        // Root of Broombridge data structure
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

            [YamlMember(Alias = "problem_description", ApplyNamingConventions = false)]
            public List<ProblemDescription> ProblemDescriptions { get; set; }

        }

        public struct ProblemDescription
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
            public List<State> InitialStates { get; set; }
        }



        public struct State
        {
            [YamlMember(Alias = "label", ApplyNamingConventions = false)]
            public string Label { get; set; }

            [YamlMember(Alias = "method", ApplyNamingConventions = false)]
            public string Method { get; set; }

            [YamlMember(Alias = "energy", ApplyNamingConventions = false)]
            public SimpleQuantity Energy { get; set; }

            /// <summary>
            /// Sparse multi-configurational data
            /// </summary>
            [YamlMember(Alias = "superposition", ApplyNamingConventions = false)]
            public List<List<string>> Superposition { get; set; }
                
            /// <summary>
            /// Coupled-cluster operator
            /// </summary>
            [YamlMember(Alias = "cluster_operator", ApplyNamingConventions = false)]
            public ClusterOperator ClusterOperator { get; set; }
        }

        /// <summary>
        /// Coupled-cluster operator
        /// </summary>
        public struct ClusterOperator
        {
            [YamlMember(Alias = "reference_state", ApplyNamingConventions = false)]
            public List<string> Reference { get; set; }

            [YamlMember(Alias = "one_body_amplitudes", ApplyNamingConventions = false)]
            public List<List<string>> OneBodyAmplitudes { get; set; }

            [YamlMember(Alias = "two_body_amplitudes", ApplyNamingConventions = false)]
            public List<List<string>> TwoBodyAmplitudes { get; set; }
        }

        /// <summary>
        /// Builds Hamiltonian from Broombridge if data is available.
        /// </summary>
        internal static OrbitalIntegralHamiltonian ToOrbitalIntegralHamiltonian(ProblemDescription broombridge)
        {
            // Add the identity terms
            var identityterm = broombridge.CoulombRepulsion.Value + broombridge.EnergyOffset.Value;
            var hamiltonian = V0_1.ToOrbitalIntegralHamiltonian(broombridge.Hamiltonian);
            //hamiltonian.Add(new OrbitalIntegral(), identityterm);
            return hamiltonian;
        }
        
        /// <summary>
        /// Parses the method field to determine the initial state.
        /// </summary>
        /// <param name="state">String in method field.</param>
        /// <returns>Initial state enum type.</returns>
        internal static StateType ParseInitialStateMethod(string state) =>
            StateTypeDictionary.ContainsKey(state) ? StateTypeDictionary[state] : StateType.Default;
        

        // Parses MCF wavefunction.
        internal static SparseMultiCFWavefunction<SpinOrbital> ToSparseMultiCFWavefunction(List<List<string>> superposition)
        {
            var outputState = new SparseMultiCFWavefunction<SpinOrbital>();
            // Todo: modify broombridge to have explicit reference state.
            foreach (var element in superposition)
            {
                // First item is the amplitude.
                double amplitude = double.Parse(element.First().ToString(), System.Globalization.CultureInfo.InvariantCulture);

                // Second to Second-last are fermion operators.
                IEnumerable<LadderOperator<SpinOrbital>> operators = element
                    .Take(element.Count() - 1)
                    .Skip(1)
                    .Select(o => V0_1.ParsePolishNotation(o.ToString()));

                var ladderSequence = new LadderSequence<SpinOrbital>(operators);

                // Sort operators to index order.
                IEnumerable<IndexOrderedSequence<SpinOrbital>> sortedOperators = ladderSequence.ToIndexOrder();
                
                // Only select the shortest term with no repeated indices.
                IndexOrderedSequence<SpinOrbital> sortedOperator = sortedOperators.ToList().OrderBy(o => o.UniqueIndices()).First();

                outputState.Set(new Complex(amplitude, 0.0), sortedOperator);
            }
            return outputState;
        }

        // Parses UCC wavefunction.
        internal static UnitaryCCWavefunction<SpinOrbital> ToUnitaryCCWavefunction(ClusterOperator clusterOperator)
        {
            var outputState = new UnitaryCCWavefunction<SpinOrbital>();

            foreach (var element in clusterOperator.OneBodyAmplitudes
                .Concat(clusterOperator.TwoBodyAmplitudes))
            {
                // First item is the amplitude.
                double amplitude = double.Parse(element.First().ToString(), System.Globalization.CultureInfo.InvariantCulture);

                // Second to last are fermion operators.
                IEnumerable<LadderOperator<SpinOrbital>> operators = element
                    .Skip(1)
                    .Select(o => V0_1.ParsePolishNotation(o.ToString()));

                var ladderSequence = new LadderSequence<SpinOrbital>(operators);

                // Terms are assumed to be in normal order.
                var sortedOperator = new IndexOrderedSequence<SpinOrbital>(operators);

                outputState.Set(new Complex(amplitude, 0.0), sortedOperator);
            }

            var reference = new List<List<string>>() { clusterOperator.Reference };

            outputState.Reference = new SingleCFWavefunction<SpinOrbital>(
               ToSparseMultiCFWavefunction(reference).Excitations.Single().Key.ToIndices()
            );

            return outputState;
        }

        // Parses all wavefunctions.
        internal static (StateType, double, object) ToWavefunction(V0_2.State state)
        {
            StateType method = ParseInitialStateMethod(state.Method);
            double energy = state.Energy.Value;
            List<List<string>> superposition = state.Superposition;
            ClusterOperator clusterOperator = state.ClusterOperator;
            object outputState;

            if (method == StateType.SparseMultiConfigurational)
            {
                outputState = ToSparseMultiCFWavefunction(superposition);
            }
            else if (method == StateType.UnitaryCoupledCluster)
            {
                outputState = ToUnitaryCCWavefunction(clusterOperator);
            }
            else
            {
                throw new System.ArgumentException($"initial state `{state.Label}` is not recognized or implemented.");
            }
            return (method, energy, outputState);
        }

    }
    #endregion
    
}
