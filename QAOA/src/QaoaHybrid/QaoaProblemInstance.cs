// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Qaoa.QaoaHybrid
{
    using System;

    /// <summary>
    /// This class is used for storing the encoding of a combinatorial optimization problem into a Hamiltonian. Currently, a Hamiltonian with locality of up to 2 is supported.
    /// </summary>
    public class QaoaProblemInstance
    {
        public double[] OneLocalHamiltonianCoefficients { get; }

        public double[] TwoLocalHamiltonianCoefficients { get; }

        public int ProblemSizeInBits { get; }

        /// <summary>
        /// Initializes a new instance of the <see cref="QaoaProblemInstance"/> class.
        /// </summary>
        /// <param name="oneLocalHamiltonianCoefficients">
        /// Given a problem encoding into a Hamiltonian, this field corresponds to coefficients of 1-local terms. Assuming that a solution to a combinatorial optimization problem can be encoded into n bits (which then corresponds to an encoding into n qubits), this array must be of length n. The i-th coefficient in an array corresponds to the coefficient of a term \sigma_i^z.
        /// </param>
        /// <param name="twoLocalHamiltonianCoefficients">
        /// Given a problem encoding into a Hamiltonian, this field corresponds to coefficients of 2-local terms. Assuming that a solution to a combinatorial optimization problem can be encoded into n bits (which then corresponds to an encoding into n qubits), this array must be of length n^2. The (i*n+j)-th coefficient in an array corresponds to the coefficient of a term \sigma_i^z\sigma_j^z (other coefficients in an array can take any value of type double).
        /// </param>
        public QaoaProblemInstance(double[] oneLocalHamiltonianCoefficients, double[] twoLocalHamiltonianCoefficients)
        {
            this.OneLocalHamiltonianCoefficients = oneLocalHamiltonianCoefficients;
            this.TwoLocalHamiltonianCoefficients = twoLocalHamiltonianCoefficients;
            this.ProblemSizeInBits = oneLocalHamiltonianCoefficients.Length;
        }

        /// <summary>
        /// Calculates the value of the objective function Hamiltonian for a binary string provided.
        /// </summary>
        /// <param name="result">
        /// A binary string. In this context it is a result that we get after measuring the QAOA state.
        /// </param>
        /// <returns>
        /// The value of the objective function Hamiltonian.
        /// </returns>
        /// <remarks>
        /// In the binary string, 0 is mapped to 1 and 1 is mapped to -1 since (-1,1) are eigenvalues of the Z operator which is currently supported in this implementation.
        /// </remarks>
        public double EvaluateHamiltonian(bool[] result)
        {
            double hamiltonianValue = 0;
            for (var i = 0; i < this.ProblemSizeInBits; i++)
            {
                hamiltonianValue += this.OneLocalHamiltonianCoefficients[i] * (1 - (2 * Convert.ToInt32(result[i])));
            }

            for (var i = 0; i < this.ProblemSizeInBits; i++)
            {
                for (var j = i + 1; j < this.ProblemSizeInBits; j++)
                {
                    hamiltonianValue += this.TwoLocalHamiltonianCoefficients[(i * this.ProblemSizeInBits) + j] * (1 - (2 * Convert.ToInt32(result[i]))) * (1 - (2 * Convert.ToInt32(result[j])));
                }
            }

            return hamiltonianValue;
        }
    }
}
