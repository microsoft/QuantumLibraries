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
    public class Hydrogen
    {
        public const double TrotterStepSize = 1.0;
        public const int TrotterOrder = 1;
        public const string GroundState = "|G>";
        // This energy is only approximately correct.
        public const double GroundStateEnergy = -1.13273749;


        public class Version_v0_1
        {
            static string filename = "Molecules/MolecularHydrogen/hydrogen_0.1.yaml";

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
        }

            public class Version_v0_2
        {
            private string TestName([CallerMemberName] string callerName = "") =>
                $"{this.GetType().FullName}.{callerName}";
            static string filename = "Molecules/MolecularHydrogen/hydrogen_0.2.yaml";

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
            
            [Fact]
            public void HighPrecisionEnergy()
            {
                var configuration = Config.Default();
                var qSharpData = Load(GroundState, configuration);
                var estEnergy = Helper.SetUpSimulation(filename, configuration, qSharpData, TrotterStepSize, TrotterOrder, 9, TestName());
                var error = GroundStateEnergy - estEnergy;

                Assert.True(Math.Abs(error) < 1e-2);
            }

            

            [Fact]
            public void EnergyUpDownIndexConvention()
            {
                var configuration = Config.Default();
                configuration.UseIndexConvention = IndexConvention.UpDown;
                var qSharpData = Load(GroundState, configuration);
                var estEnergy = Helper.SetUpSimulation(filename, configuration, qSharpData, TrotterStepSize, TrotterOrder, 9, TestName());
                var error = GroundStateEnergy - estEnergy;

                Assert.True(Math.Abs(error) < 1e-2);
            }

            [Fact]
            public void RunUnitaryCoupledCluster()
            {
                var configuration = Config.Default();
                configuration.UseIndexConvention = IndexConvention.UpDown;

                // This is a ranodm UCCSD state, not the actual one for LiH.
                var qSharpData = Load("UCCSD |G>", configuration);
                var estEnergy = Helper.SetUpSimulation(filename, configuration, qSharpData, TrotterStepSize, TrotterOrder, 9, TestName());

                var error = GroundStateEnergy - estEnergy;

                Assert.True(Math.Abs(error) < 1e-2);
            }

            [Fact]
            public void EstimateEnergyUCCSD()
            {
                var configuration = Config.Default();
                configuration.UseIndexConvention = IndexConvention.UpDown;

                // Loads a given UCCSD state
                var qSharpData = Load("UCCSD |G>", configuration);

		using (var qsim = new QuantumSimulator())
		{
		    // Estimate the energy of the molecule with UCCSD
		    var nSamples = 1000000000000000000;
                    var estEnergy = EstimateEnergy.Run(qsim, qSharpData, nSamples).Result;

		    // Compare to reference value
                    Assert.True(Math.Abs(-1.13727 - estEnergy) < 1e-3);
		}
            }

        }

    }
}
