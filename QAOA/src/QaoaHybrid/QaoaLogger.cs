namespace QAOA.QaoaHybrid
{
    using System;
    using System.IO;
    using Microsoft.Quantum.Simulation.Core;

    class QaoaLogger
    {
        private StreamWriter logger;

        public QaoaLogger(string filePath = "")
        {
            this.logger = new StreamWriter(@filePath+"hybrid_qaoa_log_" + DateTime.Now.ToString("yyyy-dd-M--HH-mm-ss") + ".txt", true);
        }

        /// # Summary
        /// Writes current values of the best fidelity and the best solution vector to a file.
        ///
        /// # Input
        /// ## bestHamiltonian
        /// Best value of a Hamiltonian so far.
        /// ## bestVector
        /// Best solution vector that generates the above value of a Hamiltonian  so far.
        public void LogCurrentBestSolution(QArray<Double> beta, QArray<Double> gamma, double bestHamiltonian, String bestVector)
        {

                this.logger.WriteLine("Current beta vector:");
                this.logger.WriteLine(beta);
                this.logger.WriteLine("Current gamma vector:");
                this.logger.WriteLine(gamma);
                this.logger.WriteLine("Current best fidelity");
                this.logger.WriteLine(bestHamiltonian);
                this.logger.WriteLine("Current best string");
                this.logger.WriteLine(bestVector);
        }

        /// # Summary
        /// Writes to a file whether an optimization finished successfully.
        ///
        /// # Input
        /// ## success
        /// A flag that indiciates whether an optimization finished successfully.
        public void LogSuccess(bool success)
        {
                this.logger.WriteLine("Was optimization successful?");
                this.logger.WriteLine(success);
                this.logger.WriteLine("##################################");

        }

        /// # Summary
        /// Closes a logger.
        public void Close()
        {
            this.logger.Close();
        }
    }
}
