// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.



using System;
using System.Runtime.CompilerServices;

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
using Microsoft.Quantum.Chemistry.JordanWigner.VQE;
using SystemTests;
using Xunit;


namespace SystemTests.Molecules
{
    public class LithiumHydride
    {

        public const double TrotterStepSize = 0.5;
        public const int TrotterOrder = 1;
        public const string GroundState = "|G>";
        public const double GroundStateEnergy = -7.881844675840115;

        public class Version_v0_1
        {
            static string filename = "Molecules/LithiumHydride/LiH_0.1.yaml";
            private string TestName([CallerMemberName] string callerName = "") =>
                $"{this.GetType().FullName}.{callerName}";

            public JordanWignerEncodingData Load(string stateName, Config config)
            {
                return Helper.GetQSharpData(filename, stateName, config);
            }

            [Fact]
            // Test classical computing Stack.
            public void LoadTest()
            {
                Load(GroundState, Config.Default());
            }

            [Fact(Skip ="Takes 2 minutes")]
            public void HighPrecisionEnergy()
            {
                var configuration = Config.Default();
                var qSharpData = Load(GroundState, configuration);
                var estEnergy = Helper.SetUpSimulation(filename, configuration, qSharpData, TrotterStepSize, TrotterOrder, 9, TestName());
                var error = GroundStateEnergy - estEnergy;

                Assert.True(Math.Abs(error) < 1e-2, "This test is probabilistic.");
            }


            [Fact]
            public void Energy()
            {
                var configuration = Config.Default();
                var qSharpData = Load(GroundState, configuration);
                var estEnergy = Helper.SetUpSimulation(filename, configuration, qSharpData, TrotterStepSize, TrotterOrder, 6, TestName());
                var error = GroundStateEnergy - estEnergy;

                Assert.True(Math.Abs(error) < 2e-1, "This test is probabilistic.");
            }

            [Fact]
            public void EnergyUpDownIndexConvention()
            {
                var configuration = Config.Default();
                configuration.UseIndexConvention = IndexConvention.UpDown;
                var qSharpData = Load(GroundState, configuration);
                var estEnergy = Helper.SetUpSimulation(filename, configuration, qSharpData, TrotterStepSize, TrotterOrder, 5);
                var error = GroundStateEnergy - estEnergy;

                Assert.True(Math.Abs(error) < 4e-1, "This test is probabilistic.");
            }

            
        }


        public class Version_v0_2
        {
            static string filename = "Molecules/LithiumHydride/LiH_0.2.yaml";
            private string TestName([CallerMemberName] string callerName = "") =>
                $"{this.GetType().FullName}.{callerName}";

            public JordanWignerEncodingData Load(string stateName, Config config)
            {
                return Helper.GetQSharpData(filename, stateName, config);
            }

            [Fact]
            public void LoadTest()
            {
                Load(GroundState, Config.Default());
            }

            [Fact(Skip = "Takes 2 minutes")]
            public void HighPrecisionEnergy()
            {
                var configuration = Config.Default();
                var qSharpData = Load(GroundState, configuration);
                var estEnergy = Helper.SetUpSimulation(filename, configuration, qSharpData, TrotterStepSize, TrotterOrder, 9, TestName());
                var error = GroundStateEnergy - estEnergy;

                Assert.True(Math.Abs(error) < 1e-2);
            }


            [Fact]
            public void Energy()
            {
                var configuration = Config.Default();
                var qSharpData = Load(GroundState, configuration);
                var estEnergy = Helper.SetUpSimulation(filename, configuration, qSharpData, TrotterStepSize, TrotterOrder, 5, TestName());
                var error = GroundStateEnergy - estEnergy;

                Assert.True(Math.Abs(error) < 4e-1);
            }

            [Fact]
            public void EnergyUpDownIndexConvention()
            {
                var configuration = Config.Default();
                configuration.UseIndexConvention = IndexConvention.UpDown;
                var qSharpData = Load(GroundState, configuration);
                var estEnergy = Helper.SetUpSimulation(filename, configuration, qSharpData, TrotterStepSize, TrotterOrder, 6, TestName());
                var error = GroundStateEnergy - estEnergy;

                Assert.True(Math.Abs(error) < 2e-1);
            }

            [Fact]
            public void RunUnitaryCoupledCluster()
            {
                var configuration = Config.Default();
                configuration.UseIndexConvention = IndexConvention.UpDown;

                // This is a random UCCSD state, not the actual one for LiH.
                var qSharpData = Load("UCCSD |E1>", configuration);
                var estEnergy = Helper.SetUpSimulation(filename, configuration, qSharpData, TrotterStepSize, TrotterOrder, 1, TestName());
            }

            [Fact]
            public void EstimateEnergyUCCSD()
            {
                var configuration = Config.Default();
                configuration.UseIndexConvention = IndexConvention.UpDown;

		// Loads a given UCCSD state
 		var qSharpData = Load("UCCSD test0", configuration);

                using (var qsim = new QuantumSimulator())
                {
                    // Estimate the energy of the molecule with UCCSD
                    var nSamples = 1000000000000000000;
                    var estEnergy = EstimateEnergy.Run(qsim, qSharpData, nSamples).Result;

                    // Compare to reference value
		    Assert.Equal(-7.8602, estEnergy, 2);
                }
            }
	}
    }
}