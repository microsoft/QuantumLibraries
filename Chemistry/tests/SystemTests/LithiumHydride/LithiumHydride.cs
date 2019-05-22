// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using System;
using System.Diagnostics;
using System.Linq;
using System.Runtime.CompilerServices;
using System.Security.Cryptography;
using System.Text;

using Microsoft.Quantum.Chemistry;
using Microsoft.Quantum.Chemistry.Broombridge;
using Microsoft.Quantum.Chemistry.Fermion;
using Microsoft.Quantum.Chemistry.JordanWigner;
using Microsoft.Quantum.Chemistry.OrbitalIntegrals;
using Microsoft.Quantum.Chemistry.QSharpFormat;
using Microsoft.Quantum.Simulation.Simulators;

using Xunit;

namespace SystemTestsLiH
{
    public class LithiumHydride
    {
        private static readonly SHA256Managed hashMethod = new SHA256Managed();

        public static JordanWignerEncodingData TestStack(string filename, Config configuration)
        {
            var broombridge = Deserializers.DeserializeBroombridge(filename).ProblemDescriptions.First();

            var orbHam = broombridge.ToOrbitalIntegralHamiltonian();

            var ferHam = orbHam.ToFermionHamiltonian(configuration.UseIndexConvention);

            var pauHam = ferHam.ToPauliHamiltonian();

            var hamiltonian = pauHam.ToQSharpFormat();

            var wavefunction = broombridge
                .ToWavefunctions(configuration.UseIndexConvention)["|G>"]
                .ToQSharpFormat();

            var qSharpData = Microsoft.Quantum.Chemistry.QSharpFormat.Convert.ToQSharpFormat(hamiltonian, wavefunction);
            return qSharpData;
        }

        public static Double SetUpLiHSimulation(string testName, string filename, Config configuration, int bits, string wavefunction = "|G>")
        {
            var qSharpData = TestStack(filename, configuration);

            // We specify the bits of precision desired in the phase estimation 
            // algorithm

            // We specify the step-size of the simulated time-evolution
            var trotterStep = 0.5;

            // Choose the Trotter integrator order
            Int64 trotterOrder = 1;

            using (var qsim = new QuantumSimulator(randomNumberGeneratorSeed: GenerateSeed(testName)))
            {

                // EstimateEnergyByTrotterization
                var (phaseEst, energyEst) = GetEnergyByTrotterization.Run(qsim, qSharpData, bits, trotterStep, trotterOrder).Result;
                var errorResult = -7.881844675840115 - energyEst;

                return errorResult;
            }
        }

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



        public class Version_v0_1
        {
            static string filename = "LithiumHydride/LiH_0.1.yaml";

            private string TestName([CallerMemberName] string callerName = "") =>
                $"{this.GetType().FullName}.{callerName}";

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

                var error = SetUpLiHSimulation(TestName(), filename, configuration, 9);

                Assert.True(Math.Abs(error) < 1e-2, "This test is probabilistic.");
            }


            [Fact]
            public void Energy()
            {
                var configuration = Config.Default();

                var error = SetUpLiHSimulation(TestName(), filename, configuration, 6);

                Assert.True(Math.Abs(error) < 1e-1, "This test is probabilistic.");
            }

            [Fact]
            public void EnergyUpDownIndexConvention()
            {
                var configuration = Config.Default();
                configuration.UseIndexConvention = SpinOrbital.IndexConvention.UpDown;

                var error = SetUpLiHSimulation(TestName(), filename, configuration, 6);

                Assert.True(Math.Abs(error) < 1e-1, "This test is probabilistic.");
            }

            
        }


        public class Version_v0_2
        {
            static string filename = "LithiumHydride/LiH_0.2.yaml";

            private string TestName([CallerMemberName] string callerName = "") =>
                $"{this.GetType().FullName}.{callerName}";

            [Fact]
            public void Load()
            {
                TestStack(filename, Config.Default());
            }

            [Fact(Skip = "Takes 2 minutes")]
            public void HighPrecisionEnergy()
            {
                var configuration = Config.Default();

                var error = SetUpLiHSimulation(TestName(), filename, configuration, 9);

                Assert.True(Math.Abs(error) < 1e-2);
            }


            [Fact]
            public void Energy()
            {
                var configuration = Config.Default();

                var error = SetUpLiHSimulation(TestName(), filename, configuration, 6);

                Assert.True(Math.Abs(error) < 1e-1);
            }

            [Fact]
            public void EnergyUpDownIndexConvention()
            {
                var configuration = Config.Default();
                configuration.UseIndexConvention = SpinOrbital.IndexConvention.UpDown;

                var error = SetUpLiHSimulation(TestName(), filename, configuration, 6);

                Assert.True(Math.Abs(error) < 1e-1);
            }

            [Fact]
            public void RunUnitaryCoupledCluster()
            {
                var configuration = Config.Default();
                configuration.UseIndexConvention = SpinOrbital.IndexConvention.UpDown;

                // This is a ranodm UCCSD state, not the actual one for LiH.
                var error = SetUpLiHSimulation(TestName(), filename, configuration, 1, "UCCSD |E1>");
            }
        }
    }
}
