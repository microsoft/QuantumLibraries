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
    class Program
    {
        static void Main(string[] args)
        {

            #region Basic Definitions

            // We start by loading the simulator that we will use to run our Q# operations.
            var qsim = new QuantumSimulator();

            // Each site of the Ising model is simulated using a single qubit. 
            var nSites = 9;

            // This coefficient is the initial coupling to the transverse field.
            var hXInitial = 1.0;

            // This coefficient is the final coupling to the transverse field.
            var hXFinal = 0.0;

            // This coefficient is the final coupling between sites.
            var jFinal = 1.0;

            // This is how long we take to sweep over the schedule parameter.
            var adiabaticTime = 10.0;

            // As we are using a Trotter-Suzuki decomposition as our simulation algorithm,
            // we will need to pick a timestep for the simulation, and the order of the 
            // integrator. The optimal timestep needs to be determined empirically, and
            // we find that the following choice works well enough.
            var trotterStepSize = 0.1;
            var trotterOrder = 2;

            // The phase estimation algorithm requires us to choose the duration of time-
            // evolution in the oracle it calls, and the bits of precision to which we
            // estimate the phase. Note that the error of the energy estimate is typically
            // rescaled by 1 / `qpeStepSize`
            var qpeStepSize = 0.1;
            var nBitsPrecision = 5;


            #endregion

            #region Ising model simulations

            // For diagnostic purposes, before we proceed to the next step, we'll print
            // out a description of the parameters we just defined.
            Console.WriteLine("\nIsing model parameters:");
            Console.WriteLine(
                $"\t{nSites} sites\n" +
                $"\t{hXInitial} initial transverse field coefficient\n" +
                $"\t{hXFinal} final transverse field coefficient\n" +
                $"\t{jFinal} final two-site coupling coefficient\n" +
                $"\t{adiabaticTime} time-interval of interpolation\n" +
                $"\t{trotterStepSize} simulation time step \n" +
                $"\t{trotterOrder} order of integrator \n" +
                $"\t{qpeStepSize} phase estimation oracle simulation time step \n" +
                $"\t{nBitsPrecision} phase estimation bits of precision\n");

            // Let us now prepare an approximate ground state of the Ising model and estimate its
            // ground state energy. This procedure is probabilistic as the quantum state
            // obtained at the end of adiabatic evolution has overlap with the ground state that
            // is less that one. Thus we repeat several times. In this case, we also know that the
            // ground state has all spins pointing in the same direction, so we print the results of
            // measuring each site after perform phase estimation to check if we were close.

            Console.WriteLine("Adiabatic state preparation of the Ising model with uniform couplings followed by phase estimation and then measurement of sites.");
            
            // Theoretical prediction of ground state energy when hXFinal is 0. 
            var energyTheory = - jFinal * ((Double)(nSites - 1));
            
            for (int rep = 0; rep < 10; rep++)
            {
                var data = IsingEstimateEnergy.Run(qsim, nSites, hXInitial, hXFinal, jFinal, adiabaticTime, trotterStepSize, trotterOrder, qpeStepSize, nBitsPrecision).Result;
                var energyEst = data.Item1;
                var stateMeasured = data.Item2.ToArray();
                Console.Write($"State: {string.Join(", ", stateMeasured.Select(x => x.ToString()).ToArray())} Energy estimate: {energyEst} vs Theory: {energyTheory}. \n");
            }

            Console.WriteLine("\n Same procedure, but using the built-in function.");
            for (int rep = 0; rep < 10; rep++)
            {
                var data = IsingEstimateEnergy_Builtin.Run(qsim, nSites, hXInitial, hXFinal, jFinal, adiabaticTime, trotterStepSize, trotterOrder, qpeStepSize, nBitsPrecision).Result;
                var phaseEst = data;
                Console.WriteLine(" Energy estimate: {0} vs Theory: {1}. ", phaseEst, energyTheory);
            }


            Console.ReadLine();
            #endregion
        }
    }
}
