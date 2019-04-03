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
        public void Energy()
        {
            var hamiltonian = FermionHamiltonian.LoadFromBroombridge(filename).First();
            var jordanWignerEncoding = JordanWignerEncoding.Create(hamiltonian);
            var qSharpData = jordanWignerEncoding.QSharpData("Greedy");


            // We specify the bits of precision desired in the phase estimation 
            // algorithm
            var bits = 8;

            // We specify the step-size of the simulated time-evolution
            var trotterStep = 0.4;

            // Choose the Trotter integrator order
            Int64 trotterOrder = 1;

            using (var qsim = new QuantumSimulator())
            {
                for (int i = 0; i < 5; i++)
                {
                    // EstimateEnergyByTrotterization
                    var (phaseEst, energyEst) = GetEnergyByTrotterization.Run(qsim, qSharpData, bits, trotterStep, trotterOrder).Result;

                    Console.WriteLine($"Rep #{i + 1}/5: Energy estimate: {energyEst}; Phase estimate: {phaseEst}");
                }
            }
        }
        

    }
}
