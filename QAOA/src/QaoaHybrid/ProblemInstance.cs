namespace Microsoft.Quantum.QAOA.QaoaHybrid
{
    using System;

    public class ProblemInstance
    {
        public double[] oneLocalHamiltonianCoefficients { get; }
        public double[] twoLocalHamiltonianCoefficients { get; }
        public int problemSizeInBits { get; }

        public ProblemInstance(double[] oneLocalHamiltonianCoefficients, double[] twoLocalHamiltonianCoefficients)
        {
            this.oneLocalHamiltonianCoefficients = oneLocalHamiltonianCoefficients;
            this.twoLocalHamiltonianCoefficients = twoLocalHamiltonianCoefficients;
            this.problemSizeInBits = oneLocalHamiltonianCoefficients.Length;
        }
    }
}
