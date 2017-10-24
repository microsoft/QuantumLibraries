using Microsoft.Quantum.Simulation.Core;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;
using System.Runtime.Serialization;
using System.Text;
using System.Threading.Tasks;
using Xunit.Sdk;

namespace Microsoft.Quantum.Canon.Tests
{
    class RandomDoubleDataAttribute : DataAttribute
    {
        private int NSamples;
        private Random RandomGenerator;
        private Double Min, Max;

        public RandomDoubleDataAttribute(int nSamples = 20, int seed = 0, Double min = 0, Double max = 1)
        {
            NSamples = nSamples;
            RandomGenerator = new Random(seed); 
            Min = min;
            Max = max;
        }
        
        public override IEnumerable<object[]> GetData(MethodInfo testMethod)
        {
            foreach (var idxSample in Enumerable.Range(0, NSamples))
            {
                yield return new[] { Min + (Max - Min) * RandomGenerator.NextDouble() as object }; 
            }
        }
    }

}
