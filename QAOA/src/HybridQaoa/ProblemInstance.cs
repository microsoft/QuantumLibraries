using System;

namespace Quantum.QAOA
{
    public class ProblemInstance
    {
        public Double[] oneLocalHamiltonianCoefficients { get; }
        public Double[] twoLocalHamiltonianCoefficients { get; }
        public int problemSizeInBits { get; }

        public ProblemInstance(Double[] oneLocalHamiltonianCoefficients, Double[] twoLocalHamiltonianCoefficients)
        {
            this.oneLocalHamiltonianCoefficients = oneLocalHamiltonianCoefficients;
            this.twoLocalHamiltonianCoefficients = twoLocalHamiltonianCoefficients;
            this.problemSizeInBits = oneLocalHamiltonianCoefficients.Length;
        }
    }
}
