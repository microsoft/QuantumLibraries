// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.



using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using Microsoft.Quantum.Chemistry;
using Microsoft.Quantum.Simulation.Core;
using Microsoft.Quantum.Simulation.Simulators;
using Microsoft.Quantum.Chemistry.Broombridge;
using Microsoft.Quantum.Chemistry.OrbitalIntegrals;
using Microsoft.Quantum.Chemistry.Fermion;
using Microsoft.Quantum.Chemistry.Pauli;
using Microsoft.Quantum.Chemistry.QSharpFormat;
using Microsoft.Quantum.Chemistry.Generic;
using Microsoft.Quantum.Chemistry.JordanWigner;
using SystemTests;
using Xunit;



namespace SystemTests.Molecules
{
    public static class Helper
    {
        public static JordanWignerEncodingData GetQSharpData(string filename, string wavefunctionLabel, Config configuration)
        {
            var broombridge = Deserializers.DeserializeBroombridge(filename).ProblemDescriptions.First();

            var orbHam = broombridge.OrbitalIntegralHamiltonian;

            var ferHam = orbHam.ToFermionHamiltonian(configuration.UseIndexConvention);

            var pauHam = ferHam.ToPauliHamiltonian();

            var hamiltonian = pauHam.ToQSharpFormat();

            Dictionary<string, FermionWavefunction<SpinOrbital>> inputStates = broombridge.Wavefunctions ?? new Dictionary<string, FermionWavefunction<SpinOrbital>>();

            // If no states are provided, use the Hartree--Fock state.
            // As fermion operators the fermion Hamiltonian are already indexed by, we now apply the desired
            // spin-orbital -> integer indexing convention.
            FermionWavefunction<int> wavefunction = inputStates.ContainsKey(wavefunctionLabel)
                ? inputStates[wavefunctionLabel].ToIndexing(configuration.UseIndexConvention)
                : ferHam.CreateHartreeFockState(broombridge.NElectrons);

            var qSharpData = Microsoft.Quantum.Chemistry.QSharpFormat.Convert.ToQSharpFormat(hamiltonian, wavefunction.ToQSharpFormat());
            return qSharpData;
        }

        public static Double SetUpSimulation(
            string filename, 
            Config configuration, 
            JordanWignerEncodingData qSharpData,
            double trotterStep,
            int trotterOrder,
            int bits
            )
        {
            
            // We specify the bits of precision desired in the phase estimation 
            // algorithm

            // We specify the step-size of the simulated time-evolution
            
            // Choose the Trotter integrator order
            
            using (var qsim = new QuantumSimulator())
            {

                // EstimateEnergyByTrotterization
                var (phaseEst, energyEst) = GetEnergyByTrotterization.Run(qsim, qSharpData, bits, trotterStep, trotterOrder).Result;
                return energyEst;
            }
        }
        

    }
}
