// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using Microsoft.Quantum.Simulation.Core;
using Microsoft.Quantum.Simulation.Simulators;

namespace Microsoft.Quantum.Samples.Ising
{
    class AdiabaticIsing
    {
        static void Main(string[] args)
        {
            var qsim = new QuantumSimulator();
            var nSites = 12;
            var trotterOrder = 2;
            var hxCoeff = 1.0;
            var jCCoeff = 1.0;


            Console.WriteLine("Ising model with {0} sites and uniform couplings", nSites);
            Console.WriteLine("Let us consider the results of fast non-adiabatic evolution from the transverse Hamiltonian to the coupling Hamiltonian.");
            var adiabaticTime = 0.1;
            var trotterStepSize = 0.01;
            for (int rep = 0; rep < 10; rep++)
            {

                
                var data = Ising1DAdiabaticAndMeasure.Run(qsim, nSites, hxCoeff, jCCoeff, adiabaticTime, trotterStepSize, trotterOrder).Result;
                var stateMeasured = data.ToArray();
                Console.Write("State: ");
                for (int idx = 0; idx < nSites; idx++)
                {
                    Console.Write("{0}  ", stateMeasured[idx]);
                }
                Console.WriteLine("");
            }

            Console.WriteLine("Observe that the zeros and ones occur randomly.");
            Console.WriteLine("As we are not in the ground state, phase estimation returns the energy of a randomly chosen eigenstate.");


            Console.WriteLine("Let us now slow down the evolution.");
            adiabaticTime = 10.0;
            trotterStepSize = 0.1;
            for (int rep = 0; rep < 10; rep++)
            {

                
                var data = Ising1DAdiabaticAndMeasure.Run(qsim, nSites, hxCoeff, jCCoeff, adiabaticTime, trotterStepSize, trotterOrder).Result;
                var stateMeasured = data.ToArray();
                Console.Write("State: ");
                for (int idx = 0; idx < nSites; idx++)
                {
                    Console.Write("{0}  ", stateMeasured[idx]);
                }
                Console.WriteLine("");
            }

            Console.WriteLine("Observe that there there now is a bias towards the sites being all equal.");
            Console.WriteLine("");
            Console.WriteLine("Studying anti-ferromagnetic coupling requires us to change the sign of jCCoeff.");
            adiabaticTime = 10.0;
            trotterStepSize = 0.1;
            jCCoeff = -1.0;
            for (int rep = 0; rep < 10; rep++)
            {


                var data = Ising1DAdiabaticAndMeasure.Run(qsim, nSites, hxCoeff, jCCoeff, adiabaticTime, trotterStepSize, trotterOrder).Result;
                var stateMeasured = data.ToArray();
                Console.Write("State: ");
                for (int idx = 0; idx < nSites; idx++)
                {
                    Console.Write("{0}  ", stateMeasured[idx]);
                }
                Console.WriteLine("");
            }


            Console.ReadLine();
        }
    }
}
