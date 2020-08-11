// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.QAOA.QaoaHybrid
{
    using System;

    public class OptimalSolution
    {

        public string OptimalVector { get; }
        public double OptimalValue { get; }
        public double[] OptimalBeta { get; }
        public double[] OptimalGamma { get; }

        public OptimalSolution(string optimalVector, double optimalValue, double[] optimalBeta, double[] optimalGamma)
        {
            this.OptimalVector = optimalVector;
            this.OptimalValue = optimalValue;
            this.OptimalBeta = optimalBeta;
            this.OptimalGamma = optimalGamma;

        }

    }
}
