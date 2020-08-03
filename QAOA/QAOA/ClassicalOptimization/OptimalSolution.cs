using System;
using System.Collections.Generic;
using System.Text;

namespace QAOA.ClassicalOptimization
{
    public class OptimalSolution
    {

        public String optimalVector { get; }
        public Double optimalValue { get; }
        public Double[] optimalBeta { get; }
        public Double[] optimalGamma { get; }

        public OptimalSolution(String optimalVector, Double optimalValue, Double[] optimalBeta, Double[] optimalGamma)
        {
            this.optimalVector = optimalVector;
            this.optimalValue = optimalValue;
            this.optimalBeta = optimalBeta;
            this.optimalGamma = optimalGamma;

        }

    }
}
