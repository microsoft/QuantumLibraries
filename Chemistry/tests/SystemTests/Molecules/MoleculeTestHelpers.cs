// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.



using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Diagnostics;
using System.Runtime.CompilerServices;
using System.Security.Cryptography;


using Microsoft.Quantum.Chemistry;
using Microsoft.Quantum.Simulation.Core;
using Microsoft.Quantum.Simulation.Simulators;
using Microsoft.Quantum.Chemistry.Broombridge;
using Microsoft.Quantum.Chemistry.OrbitalIntegrals;
using Microsoft.Quantum.Chemistry.Fermion;
using Microsoft.Quantum.Chemistry.Paulis;
using Microsoft.Quantum.Chemistry.QSharpFormat;
using Microsoft.Quantum.Chemistry.Generic;
using Microsoft.Quantum.Chemistry.JordanWigner;
using SystemTests;
using Xunit;
using System.IO;

namespace SystemTests.Molecules
{
    public static class Helper
    {

        private static readonly SHA256 hashMethod = SHA256.Create();

        /// <summary>
        /// Returns a seed to use for the test run based on the class
        /// </summary>
        public static uint? GenerateSeed(string testName)
        {
            byte[] bytes = Encoding.Unicode.GetBytes(testName);
            byte[] hash = hashMethod.ComputeHash(bytes);
            uint seed = BitConverter.ToUInt32(hash, 0);

            string msg = $"Using generated seed: (\"{testName}\",{ seed })";
            Console.WriteLine(msg);
            Debug.WriteLine(msg);

            return seed;
        }


        public static JordanWignerEncodingData GetQSharpData(string filename, string wavefunctionLabel, Config configuration)
        {
            using var reader = File.OpenText(filename);
            var broombridge = BroombridgeSerializer.Deserialize(reader).First();

            var orbHam = broombridge.OrbitalIntegralHamiltonian;

            var ferHam = orbHam.ToFermionHamiltonian(configuration.UseIndexConvention);

            var pauHam = ferHam.ToPauliHamiltonian();

            var hamiltonian = pauHam.ToQSharpFormat();

            var inputStates = broombridge.InitialStates ?? new Dictionary<string, FermionWavefunction<SpinOrbital>>();

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
            int bits,
            string testName = "Default"
            )
        {
            
            // We specify the bits of precision desired in the phase estimation 
            // algorithm

            // We specify the step-size of the simulated time-evolution
            
            // Choose the Trotter integrator order
            
            using (var qsim = new QuantumSimulator(randomNumberGeneratorSeed: GenerateSeed(testName)))
            {

                // EstimateEnergyByTrotterization
                var (phaseEst, energyEst) = GetEnergyByTrotterization.Run(qsim, qSharpData, bits, trotterStep, trotterOrder).Result;
                return energyEst;
            }
        }
        

    }
}



