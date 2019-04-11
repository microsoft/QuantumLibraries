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


using Microsoft.Quantum.Chemistry.Fermion;
using Microsoft.Quantum.Chemistry.OrbitalIntegrals;

namespace Microsoft.Quantum.Chemistry.Broombridge
{
    /// <summary>
    /// Latest Broombridge format.
    /// </summary>
    public class Data
    {
        /// <summary>
        /// Unparsed deserialized Broombridge data.
        /// </summary>
        public V0_2.Data Raw { get; set; }

        // Root of Broombridge data structure

        /// <summary>
        /// URL to schema defining this version of Broombridge.
        /// </summary>
        public string Schema { get; set; }

        /// <summary>
        /// Broombridge version number
        /// </summary>
        public VersionNumber VersionNumber { get; set; }

        /// <summary>
        /// Collection of electronic structure problem.
        /// </summary>
        public List<ProblemDescription> ProblemDescriptions { get; set; } = new List<ProblemDescription>();

        /// <summary>
        /// Electronic structure problem instance.
        /// </summary>
        public struct ProblemDescription
        {
            /// <summary>
            /// Identity term of the Hamilonian.
            /// </summary>
            public double EnergyOffset { get; set; }

            /// <summary>
            /// Number of orbitals.
            /// </summary>
            public int NOrbitals { get; set; }

            /// <summary>
            /// Number of electroncs.
            /// </summary>
            public int NElectrons { get; set; }

            /// <summary>
            /// Hamiltonian represented by orbital integrals.
            /// </summary>
            public OrbitalIntegralHamiltonian OrbitalIntegralHamiltonian { get; set; }

            /// <summary>
            /// Collection of trial wavefunctions.
            /// </summary>
            public Dictionary<string,FermionWavefunction<SpinOrbital>> Wavefunctions { get; set; }

        }

        /// <summary>
        /// Parses deserialized Broombridge
        /// </summary>
        /// <param name="broombridgeV0_2">Broombridge data structure.</param>
        internal Data(Broombridge.V0_2.Data broombridgeV0_2)
        {
            Raw = broombridgeV0_2;
            Schema = broombridgeV0_2.Schema;
            VersionNumber = VersionNumber.v0_2;

            foreach (var problem in Raw.ProblemDescriptions)
            {
                var problemDescription = new ProblemDescription();
                problemDescription.EnergyOffset = problem.EnergyOffset.Value + problem.CoulombRepulsion.Value;
                problemDescription.NElectrons = problem.NElectrons;
                problemDescription.NOrbitals = problem.NOrbitals;
                problemDescription.OrbitalIntegralHamiltonian = V0_2.ToOrbitalIntegralHamiltonian(problem);
                problemDescription.Wavefunctions = new Dictionary<string, FermionWavefunction<SpinOrbital>>();
                foreach (var initialState in problem.InitialStates)
                {
                    var state = new FermionWavefunction<SpinOrbital>();

                    var (method, energy, outputState) = V0_2.ToWavefunction(initialState);
                    
                    state.Method = V0_2.ParseInitialStateMethod(initialState.Method);
                    state.Energy = energy;

                    if(state.Method == StateType.SparseMultiConfigurational)
                    {
                        state.MCFData = (SparseMultiCFWavefunction<SpinOrbital>) outputState;
                    }
                    else if(state.Method == StateType.UnitaryCoupledCluster)
                    {
                        state.UCCData = (UnitaryCCWavefunction<SpinOrbital>) outputState;
                    }

                    problemDescription.Wavefunctions.Add(initialState.Label, state);
                }
                ProblemDescriptions.Add(problemDescription);
            }
        }
    }
}