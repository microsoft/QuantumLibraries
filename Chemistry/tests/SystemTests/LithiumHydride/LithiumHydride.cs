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

using Xunit;

namespace SystemTests
{
    using static FermionTermType.Common;
    using FermionTerm = FermionTerm;
    using FermionTermType = FermionTermType;
    public class LithiumHydride
    {
        static string filename = "LithiumHydride/lih_sto-3g_fci_1.624.yaml";

        [Fact]
        public void Load()
        {
            var hamiltonian = FermionHamiltonian.LoadFromBroombridge(filename).First();
            var jordanWignerEncoding = JordanWignerEncoding.Create(hamiltonian);
            var qSharpData = jordanWignerEncoding.QSharpData("|G>");
        }

        [Fact]
        public void HighPrecisionEnergy()
        {
            var configuration = Config.Default();

            var error = SetUpLiHSimulation(configuration, 9);

            Assert.True(Math.Abs(error) < 1e-2);
        }


        [Fact]
        public void Energy()
        {
            var configuration = Config.Default();

            var error = SetUpLiHSimulation(configuration, 6);

            Assert.True(Math.Abs(error) < 1e-1);
        }

        [Fact]
        public void EnergyUpDownIndexConvention()
        {
            var configuration = Config.Default();
            configuration.indexConvention = SpinOrbital.Config.IndexConvention.Type.UpDown;

            var error = SetUpLiHSimulation(configuration, 6);

            Assert.True(Math.Abs(error) < 1e-1);                
        }

        public Double SetUpLiHSimulation(Config configuration, int bits)
        {
            var hamiltonian = FermionHamiltonian.LoadFromBroombridge(filename, configuration).First();
            var jordanWignerEncoding = JordanWignerEncoding.Create(hamiltonian);
            var qSharpData = jordanWignerEncoding.QSharpData("|G>");

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

    }
}
