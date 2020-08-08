namespace Microsoft.Quantum.QAOA.QaoaHybrid
{
    using System;

    public class OptimalSolution
    {

        public string optimalVector { get; }
        public double optimalValue { get; }
        public double[] optimalBeta { get; }
        public double[] optimalGamma { get; }

        public OptimalSolution(string optimalVector, double optimalValue, double[] optimalBeta, double[] optimalGamma)
        {
            this.optimalVector = optimalVector;
            this.optimalValue = optimalValue;
            this.optimalBeta = optimalBeta;
            this.optimalGamma = optimalGamma;

        }

    }
}
