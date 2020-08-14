// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.QAOA.QaoaHybrid
{
    /// <summary>
    /// This class stores a solution that is found by a hybrid QAOA. It includes a boolean array that encodes the solution, a corresponding expected value of the objective function Hamiltonian and corresponding beta and gamma parameters that are input to a QAOA.
    /// </summary>
    /// <remarks>
    /// Note that given the nature of a QAOA and a hybrid QAOA, a solution produced by the algorithm is not necessarily feasible and not necessarily optimal.
    /// </remarks>
    public class Solution
    {

        public bool[] SolutionVector { get; set; }

        public double SolutionHamiltonianValue { get; set; }

        public double[] SolutionBeta { get; set; }

        public double[] SolutionGamma { get; set; }

        /// <summary>
        /// Initializes a new instance of the <see cref="Solution"/> class.
        /// </summary>
        /// <param name="solutionVector">
        /// A vector that is a solution of a combinatorial optimization problem found by a hybrid QAOA.
        /// </param>
        /// <param name="solutionHamiltonianValue">
        /// An expected value of an objective function Hamiltonian that corresponds to the solution stored in solutionVector.
        /// </param>
        /// <param name="solutionBeta">
        /// Values of beta parameters that correspond to the solution stored in solutionVector.
        /// </param>
        /// <param name="solutionGamma">
        /// Values of gamma parameters that corresponds to the solution stored in solutionVector.
        /// </param>
        public Solution(bool[] solutionVector, double solutionHamiltonianValue, double[] solutionBeta, double[] solutionGamma)
        {
            this.SolutionVector = solutionVector;
            this.SolutionHamiltonianValue = solutionHamiltonianValue;
            this.SolutionBeta = solutionBeta;
            this.SolutionGamma = solutionGamma;

        }

    }
}
