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
                var estEnergy = Helper.SetUpSimulation(filename, configuration, qSharpData, TrotterStepSize, TrotterOrder, 9);
                var error = GroundStateEnergy - estEnergy;

                Assert.True(Math.Abs(error) < 1e-2);
            }

            

            [Fact]
            public void EnergyUpDownIndexConvention()
            {
                var configuration = Config.Default();
                configuration.UseIndexConvention = IndexConvention.UpDown;
                var qSharpData = Load(GroundState, configuration);
                var estEnergy = Helper.SetUpSimulation(filename, configuration, qSharpData, TrotterStepSize, TrotterOrder, 9);
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
                var estEnergy = Helper.SetUpSimulation(filename, configuration, qSharpData, TrotterStepSize, TrotterOrder, 9);

                var error = GroundStateEnergy - estEnergy;

                Assert.True(Math.Abs(error) < 1e-2);
            }


        }

    }
}
