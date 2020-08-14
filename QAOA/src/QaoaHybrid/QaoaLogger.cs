// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

using System.Linq;

namespace Microsoft.Quantum.QAOA.QaoaHybrid
{
    using System;
    using System.IO;
    using Microsoft.Quantum.Simulation.Core;

    /// <summary>
    /// This class provides a simple capability for logging intermediate steps of a hybrid QAOA to a text file.
    /// </summary>
    class QaoaLogger
    {
        private readonly StreamWriter logger;

        public QaoaLogger()
        {

            this.logger = new StreamWriter("hybrid_qaoa_log_" + DateTime.Now.ToString("yyyy-dd-M--HH-mm-ss") + ".txt", true);
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

                this.logger.WriteLine("Current beta vector:");
                this.logger.WriteLine(beta);
                this.logger.WriteLine("Current gamma vector:");
                this.logger.WriteLine(gamma);
                this.logger.WriteLine("Current best expected value of a Hamiltonian:");
                this.logger.WriteLine(bestHamiltonian);
                this.logger.WriteLine("Current best solution vector:");
                var bestSolutionVector = string.Join(", ", bestVector.Select(x => x.ToString()));
                this.logger.WriteLine("[" + bestSolutionVector + "]");
        }

        /// <summary>
        /// Writes to a file whether an optimization finished successfully.
        /// </summary>
        /// <param name="success">
        /// A flag that indicates whether an optimization finished successfully.
        /// </param>
        public void LogSuccess(bool success)
        {
                this.logger.WriteLine("Was optimization successful?");
                this.logger.WriteLine(success);
                this.logger.WriteLine("##################################");

        }

        /// <summary>
        /// Closes a logger.
        /// </summary>
        public void Close()
        {
            this.logger.Close();
        }
    }
}
