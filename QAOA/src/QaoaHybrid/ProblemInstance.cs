// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.QAOA.QaoaHybrid
{
    using System;

    public class ProblemInstance
    {
        public double[] OneLocalHamiltonianCoefficients { get; }
        public double[] TwoLocalHamiltonianCoefficients { get; }
        public int ProblemSizeInBits { get; }

        public ProblemInstance(double[] oneLocalHamiltonianCoefficients, double[] twoLocalHamiltonianCoefficients)
        {
            this.OneLocalHamiltonianCoefficients = oneLocalHamiltonianCoefficients;
            this.TwoLocalHamiltonianCoefficients = twoLocalHamiltonianCoefficients;
            this.ProblemSizeInBits = oneLocalHamiltonianCoefficients.Length;
        }
    }
}
