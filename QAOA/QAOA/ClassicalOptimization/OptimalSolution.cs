using System;
using System.Collections.Generic;
using System.Text;

namespace QAOA.ClassicalOptimization
{
    public class OptimalSolution
    {

        public String OptimalVector { get; }
        public Double OptimalValue { get; }
        public Double[] OptimalBeta { get; }
        public Double[] OptimalGamma { get; }

        public OptimalSolution(String optimalVector, Double optimalValue, Double[] optimalBeta, Double[] optimalGamma)
        {
            OptimalVector = optimalVector;
            OptimalValue = optimalValue;
            OptimalBeta = optimalBeta;
            OptimalGamma = optimalGamma;

        }

    }
}
