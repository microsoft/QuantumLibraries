// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Qaoa.QaoaHybrid
{
    /// <summary>
    /// This class stores a solution that is found by a hybrid QAOA. It includes a boolean array that encodes the solution, a corresponding expected value of the objective function Hamiltonian and corresponding beta and gamma parameters that are input to a QAOA.
    /// </summary>
    /// <remarks>
    /// Note that given the nature of a QAOA and a hybrid QAOA, a solution produced by the algorithm is not necessarily feasible and not necessarily optimal.
    /// </remarks>
    public class QaoaSolution
    {

        public bool[] SolutionVector { get; set; }

        public double SolutionHamiltonianValue { get; set; }

        public QaoaParameters SolutionQaoaParameters { get; set; }

        /// <summary>
        /// Initializes a new instance of the <see cref="QaoaSolution"/> class.
        /// </summary>
        /// <param name="solutionVector">
        /// A vector that is a solution of a combinatorial optimization problem found by a hybrid QAOA.
        /// </param>
        /// <param name="solutionHamiltonianValue">
        /// An expected value of an objective function Hamiltonian that corresponds to the solution stored in solutionVector.
        /// </param>
        /// <param name="solutionQaoaParameters">
        /// Values of beta and gamma parameters that correspond to the solution stored in solutionVector.
        /// </param>
        public QaoaSolution(bool[] solutionVector, double solutionHamiltonianValue, QaoaParameters solutionQaoaParameters)
        {
            SolutionVector = solutionVector;
            SolutionHamiltonianValue = solutionHamiltonianValue;
            SolutionQaoaParameters = solutionQaoaParameters;

        }

    }
}
