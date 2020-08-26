// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Qaoa.QaoaHybrid
{
    using System;
    using System.IO;
    using Microsoft.Quantum.Simulation.Core;
    using System.Linq;

    /// <summary>
    /// This class provides a simple capability for logging intermediate steps of a hybrid QAOA to a text file.
    /// </summary>
    class QaoaLogger : IDisposable
    {
        private readonly StreamWriter logger;

        public QaoaLogger()
        {

            logger = new StreamWriter("hybrid_qaoa_log_" + DateTime.Now.ToString("yyyy-dd-M--HH-mm-ss") + ".txt", true);
        }

        /// <summary>
        /// Writes current values of the best fidelity and the best solution vector to a file.
        /// </summary>
        /// <param name="beta">
        /// Best beta vector so far.
        /// </param>
        /// <param name="gamma">
        /// Best gamma vector so far.
        /// </param>
        /// <param name="bestHamiltonian">
        /// Best value of a Hamiltonian so far.
        /// </param>
        /// <param name="bestVector">
        /// Best solution vector that generates the above value of a Hamiltonian  so far.
        /// </param>
        public void LogCurrentBestSolution(QArray<double> beta, QArray<double> gamma, double bestHamiltonian, bool[] bestVector)
        {

                logger.WriteLine("Current beta vector:");
                logger.WriteLine(beta);
                logger.WriteLine("Current gamma vector:");
                logger.WriteLine(gamma);
                logger.WriteLine("Current best expected value of a Hamiltonian:");
                logger.WriteLine(bestHamiltonian);
                logger.WriteLine("Current best solution vector:");
                var bestSolutionVector = string.Join(", ", bestVector.Select(x => x.ToString()));
                logger.WriteLine("[" + bestSolutionVector + "]");
        }

        /// <summary>
        /// Writes to a file whether an optimization finished successfully.
        /// </summary>
        /// <param name="success">
        /// A flag that indicates whether an optimization finished successfully.
        /// </param>
        public void LogSuccess(bool success)
        {
                logger.WriteLine("Was optimization successful?");
                logger.WriteLine(success);
                logger.WriteLine("##################################");

        }

        protected virtual void Dispose(bool disposing)
        {

            if (disposing)
            {
                logger?.Dispose();
            }
        }

        /// <summary>
        /// Closes a logger.
        /// </summary>
        public void Dispose()
        {
            Dispose(true);
        }
    }
}
