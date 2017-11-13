using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using Microsoft.Quantum.Simulation.Core;
using Microsoft.Quantum.Simulation.Simulators;

// FIXME: add prose in comments. See TeleportationSample for a comparison.

namespace Microsoft.Quantum.Samples.Ising
{
    class IsingPhaseEstimation
    {
        static void Main(string[] args)
        {
            var qsim = new QuantumSimulator();

            var nSites = 8;
            var trotterOrder = 2;
            var hXInitial = 1.0;
            var hXFinal = 0.0;
            var jZFinal = 1.0;
            var adiabaticTime = 100.0;
            var trotterStepSize = 0.1;

            var qpeStepSize = 0.1;
            var nBitsPrecision = 8;

            // Theoretical prediction of ground state energy when hXFinal is 0. 
            var energyTheory = - jZFinal * ((Double)(nSites - 1));

            Console.WriteLine("Ising model with uniform couplings.");
            Console.WriteLine("Adiabatic evolution followed by quatum phase estimation and then measurement of sites.");
            for (int rep = 0; rep < 10; rep++)
            {
                var data = IsingEstimateEnergy.Run(qsim, nSites, hXInitial, hXFinal, jZFinal, adiabaticTime, trotterStepSize, trotterOrder, qpeStepSize, nBitsPrecision).Result;
                var phaseEst = data.Item1;
                var stateMeasured = data.Item2.ToArray();
                Console.Write("State: ");
                for (int idx = 0; idx < nSites; idx++)
                {
                        Console.Write("{0} ", stateMeasured[idx]);
                }
                Console.WriteLine(" Energy estimate: {0} vs Theory: {1}. ", phaseEst, energyTheory);
            }

            Console.WriteLine("");
            Console.WriteLine("Adiabatic evolution followed by quatum phase estimation using built-in function.");
            for (int rep = 0; rep < 10; rep++)
            {
                var data = IsingEstimateEnergy_Builtin.Run(qsim, nSites, hXInitial, hXFinal, jZFinal, adiabaticTime, trotterStepSize, trotterOrder, qpeStepSize, nBitsPrecision).Result;
                var phaseEst = data;
                Console.WriteLine(" Energy estimate: {0} vs Theory: {1}. ", phaseEst, energyTheory);
            }


                Console.ReadLine();
        }
    }
}
