// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.QAOA.QaoaHybrid
{

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
