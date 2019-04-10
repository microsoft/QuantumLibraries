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

using Xunit;



namespace SystemTestsLiH
{
    public class LithiumHydride
    {
        public static JordanWignerEncodingData TestStack(string filename, Config configuration)
        {
            var broombridge = Deserializers.DeserializeBroombridge(filename).ProblemDescriptions.First();

            var orbHam = broombridge.OrbitalIntegralHamiltonian;

            var ferHam = orbHam.ToFermionHamiltonian(configuration.UseIndexConvention);

            var pauHam = ferHam.ToPauliHamiltonian();

            var hamiltonian = pauHam.ToQSharpFormat();

            var wavefunction = broombridge
                .Wavefunctions["|G>"].ToIndexing(configuration.UseIndexConvention)
                .ToQSharpFormat();

            var qSharpData = Microsoft.Quantum.Chemistry.QSharpFormat.Convert.ToQSharpFormat(hamiltonian, wavefunction);
            return qSharpData;
        }

        public static Double SetUpLiHSimulation(string filename, Config configuration, int bits, string wavefunction = "|G>")
        {
            var qSharpData = TestStack(filename, configuration);

            // We specify the bits of precision desired in the phase estimation 
            // algorithm

            // We specify the step-size of the simulated time-evolution
            var trotterStep = 0.5;

            // Choose the Trotter integrator order
            Int64 trotterOrder = 1;

            using (var qsim = new QuantumSimulator())
            {

                // EstimateEnergyByTrotterization
                var (phaseEst, energyEst) = GetEnergyByTrotterization.Run(qsim, qSharpData, bits, trotterStep, trotterOrder).Result;
                var errorResult = -7.881844675840115 - energyEst;

                return errorResult;
            }
        }


        public class Version_v0_1
        {
            static string filename = "LithiumHydride/LiH_0.1.yaml";

            [Fact]
            // Test classical computing Stack.
            public void Load()
            {

                TestStack(filename, Config.Default());
            }

            [Fact(Skip ="Takes 2 minutes")]
            public void HighPrecisionEnergy()
            {
                var configuration = Config.Default();

                var error = SetUpLiHSimulation(filename, configuration, 9);

                Assert.True(Math.Abs(error) < 1e-2, "This test is probabilistic.");
            }


            [Fact]
            public void Energy()
            {
                var configuration = Config.Default();

                var error = SetUpLiHSimulation(filename, configuration, 7);

                Assert.True(Math.Abs(error) < 1e-1, "This test is probabilistic.");
            }

            [Fact]
            public void EnergyUpDownIndexConvention()
            {
                var configuration = Config.Default();
                configuration.UseIndexConvention = IndexConvention.UpDown;

                var error = SetUpLiHSimulation(filename, configuration, 6);

                Assert.True(Math.Abs(error) < 2e-1, "This test is probabilistic.");
            }

            
        }


        public class Version_v0_2
        {
            static string filename = "LithiumHydride/LiH_0.2.yaml";

            [Fact]
            public void Load()
            {
                TestStack(filename, Config.Default());
            }

            [Fact(Skip = "Takes 2 minutes")]
            public void HighPrecisionEnergy()
            {
                var configuration = Config.Default();

                var error = SetUpLiHSimulation(filename, configuration, 9);

                Assert.True(Math.Abs(error) < 1e-2);
            }


            [Fact]
            public void Energy()
            {
                var configuration = Config.Default();

                var error = SetUpLiHSimulation(filename, configuration, 6);

                Assert.True(Math.Abs(error) < 2e-1);
            }

            [Fact]
            public void EnergyUpDownIndexConvention()
            {
                var configuration = Config.Default();
                configuration.UseIndexConvention = IndexConvention.UpDown;

                var error = SetUpLiHSimulation(filename, configuration, 7);

                Assert.True(Math.Abs(error) < 1e-1);
            }

            [Fact]
            public void RunUnitaryCoupledCluster()
            {
                var configuration = Config.Default();
                configuration.UseIndexConvention = IndexConvention.UpDown;

                // This is a ranodm UCCSD state, not the actual one for LiH.
                var error = SetUpLiHSimulation(filename, configuration, 1, "UCCSD |E1>");
            }
        }

    }
}
