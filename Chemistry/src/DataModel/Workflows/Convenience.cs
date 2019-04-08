// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.
using System;
using System.Linq;
using System.Collections.Generic;

using Microsoft.Quantum.Simulation.Core;

using Microsoft.Quantum.Chemistry.Broombridge;
using Microsoft.Quantum.Chemistry.OrbitalIntegrals;
using Microsoft.Quantum.Chemistry.Fermion;
using Microsoft.Quantum.Chemistry.Pauli;
using Microsoft.Quantum.Chemistry.QSharpFormat;
using Microsoft.Quantum.Chemistry.JordanWigner;
using Microsoft.Quantum.Chemistry.Generic;
using Microsoft.Quantum.Chemistry.LadderOperators;

namespace Microsoft.Quantum.Chemistry
{
    /// <summary>
    /// This class contains convenience functions for simulating electronic structure problems.
    /// </summary>
    public static class Convenience
    {
        public class ProblemContainer
        {
            public class Config { }
            //public readonly Config.IndexConvention.LadderType IndexConvention;

            // For now, this only has a few options.
            // Hamiltonian data structures.
            public OrbitalIntegralHamiltonian orbitalIntegralHamiltonian = new OrbitalIntegralHamiltonian();
            public FermionHamiltonian fermionHamiltonian = new FermionHamiltonian();
            public PauliHamiltonian pauliHamiltonian = new PauliHamiltonian();

            // Wavefunction data structure.
            //public Dictionary<string, InputState> wavefunctions = new Dictionary<string, InputState>();

            // QSharp data structure.
            public JordanWignerEncodingData qSharpData;
            
            // Additional data
            public int NOrbitals = 0;
            public int NElectrons = 0;
            public string MiscellaneousInformation;

            /// <summary>
            ///     A label for this particular Hamiltonian.
            ///     Can be used to identify the Hamiltonian out of set
            ///     loaded from the same file.
            /// </summary>
            public string Name { get; set; } = "<unknown>";

        }

        /// <summary>
        ///      Loads a Hamiltonian from integral data represented
        ///      in Broombridge format.
        ///      Please see the <a href="https://docs.microsoft.com/quantum/libraries/chemistry/schema/spec">
        ///      for further details about the
        ///      format parsed by this method.
        /// </summary>
        /// <param name="filename">The name of the file to be loaded.</param>
        /// <returns>
        ///      An instance of <see cref="FermionHamiltonian"/> representing the
        ///      data contained in <paramref name="filename"/>.
        /// </returns>
        public static IEnumerable<FermionHamiltonian> LoadFromBroombridge(
            string filename,
            SpinOrbital.IndexConvention indexConvention)
        {
            var broombridgeData = Deserializers.DeserializeBroombridge(filename);
            
            IEnumerable<CurrentVersion.ProblemDescription> problemData = broombridgeData.ProblemDescriptions;

            // Create electronic structure Hamiltonian
            var fermionHamiltonians = problemData
                .Select(o => o
                .CreateOrbitalIntegralHamiltonian()
                .ToFermionHamiltonian(indexConvention));

            return fermionHamiltonians;
        }

        /// <summary>
        /// This approximates the Hamiltonian ground state by a greedy algorithm  
        /// that minimizes only the PP term energies. If there are no PP terms,
        /// states will be occupied in lexicographic order.
        /// </summary>
        /// <returns>
        /// Greedy trial state for minimizing Hamiltonian diagonal one-electron energy.
        /// </returns>
        public static InputState GreedyStatePreparation(this FermionHamiltonian hamiltonian, int nElectrons)
        {
            InputState state = new InputState();
            WavefunctionFermionSCF greedyState = hamiltonian.GreedyStatePreparationSCF(nElectrons);

            state.Label = "Greedy";

            // Currently not used
            // state.reference = null;

            state.type = StateType.SparseMultiConfigurational;

            state.Superposition.Add(((1.0, 0.0), new IndexOrderedLadderSequence(greedyState.GetLadderSequence())));

            return state;
        }

        /*
        /// <summary>
        /// Extracts only the required information from a Broombridge problem instance.
        /// </summary>
        /// <param name="broombridgeProblem">A Broombridge problem description.</param>
        /// <param name="indexConvention">The indexing convention used to map a spin-orbital indicx to a single integer.</param>
        public BroombridgeTyped(Broombridge.Current.ProblemDescription broombridgeProblem, SpinOrbital.Config.IndexConvention.IndexConvention indexConvention = SpinOrbital.Config.IndexConvention.Default)
        {
            IndexConvention = indexConvention;
            NOrbitals = broombridgeProblem.NOrbitals;
            NElectrons = broombridgeProblem.NElectrons;

            IdentityTerm = broombridgeProblem.CoulombRepulsion.Value + broombridgeProblem.EnergyOffset.Value;

            InitialStates = broombridgeProblem.InitialStates.ToDictionary(
                o => o.Label,
                o => ParseInitialState(o, indexConvention)
                );
        }


        // Make LoadFromBroombridge with hamiltonian & states.
        */
    }
}