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
using Microsoft.Quantum.Chemistry.Broombridge;
using Microsoft.Quantum.Chemistry.OrbitalIntegrals;
using Microsoft.Quantum.Chemistry.Fermion;
using Microsoft.Quantum.Chemistry.Pauli;
using Microsoft.Quantum.Chemistry.Generic;

namespace Microsoft.Quantum.Chemistry.Broombridge
{
    /// <summary>
    /// Latest Broombridge format.
    /// </summary>
    public class Data
    {
        // Raw data.
        public V0_2.Data Raw { get; set; }

        // Root of Broombridge data structure

        public string Schema { get; set; }
        public VersionNumber VersionNumber { get; set; }
        public List<ProblemDescription> ProblemDescriptions = new List<ProblemDescription>();


        public struct ProblemDescription
        {
            public double EnergyOffset { get; set; }
            public int NOrbitals { get; set; }
            public int NElectrons { get; set; }
            public OrbitalIntegralHamiltonian OrbitalIntegralHamiltonian { get; set; }
            public Dictionary<string,FermionWavefunction<SpinOrbital>> Wavefunctions { get; set; }

        }
        

        

        public Data(Broombridge.V0_2.Data broombridgeV0_2)
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
                    state.Energy = initialState.Energy.Value;

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