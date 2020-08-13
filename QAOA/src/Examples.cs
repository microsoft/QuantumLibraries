// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Quantum.QAOA
{
    using Microsoft.Quantum.QAOA.QaoaHybrid;
    using System;


    class Examples
    {
        static void Main(string[] args)
        {
            //PARAMETERS
            var numberOfIterations = 70;
            var p = 3;
            var numberOfRandomStartingPoints = 1;

            //EXAMPLES

            //Quantum Santa (http://quantumalgorithmzoo.org/traveling_santa/)
            double[] dtx = { 0.619193, 0.742566, 0.060035, -1.568955, 0.045490 };
            double[] dtz = { 3.182203, -1.139045, 0.221082, 0.537753, -0.417222 };
            double[] oneLocalHamiltonianCoefficients = { 4 * 20 - 0.5 * 4.7, 4 * 20 - 0.5 * 9.09, 4 * 20 - 0.5 * 9.03, 4 * 20 - 0.5 * 5.70, 4 * 20 - 0.5 * 8.02, 4 * 20 - 0.5 * 1.71 };
            var twoLocalHamiltonianCoefficients = new[]
            {
                40.0, 40.0, 20.0, 40.0, 40.0, 40.0,
                40.0, 40.0, 40.0, 20.0, 40.0, 40.0,
                40.0, 40.0, 40.0, 40.0, 40.0, 40.0,
                40.0, 40.0, 40.0, 40.0, 40.0, 40.0,
                40.0, 40.0, 40.0, 40.0, 40.0, 20.0,
                40.0, 40.0, 40.0, 40.0, 40.0, 40.0,
            };
            var quantumSanta = new ProblemInstance(oneLocalHamiltonianCoefficients, twoLocalHamiltonianCoefficients);

            double[] segmentCosts = { 4.70, 9.09, 9.03, 5.70, 8.02, 1.71 };

            
            //MaxCut (medium.com/mdr-inc/qaoa-maxcut-using-blueqat-aaf33038f46e)
            oneLocalHamiltonianCoefficients = new Double[] { 0,0,0,0,0};
            twoLocalHamiltonianCoefficients = new Double[]{ 0,1,0,1,0,
                               0,0,1,0,0,
                               0,0,0,1,1,
                               0,0,0,0,1,
                               0,0,0,0,0};
            var maxCut1 = new ProblemInstance(oneLocalHamiltonianCoefficients, twoLocalHamiltonianCoefficients);

            //Rigetti MaxCut unit tests
            oneLocalHamiltonianCoefficients = new Double[]{-0.5,0,-1,0.5};
            twoLocalHamiltonianCoefficients = new Double[]{0,1,2,0,
                              0,0,0.5,0,
                              0,0,0,2.5,
                              0,0,0,0};
            var maxCut2 = new ProblemInstance(oneLocalHamiltonianCoefficients, twoLocalHamiltonianCoefficients);

            oneLocalHamiltonianCoefficients = new Double[] { 0.8, -0.5 };
            twoLocalHamiltonianCoefficients = new Double[]{ 0, -1, 0, 0};
            var maxCut3 = new ProblemInstance(oneLocalHamiltonianCoefficients, twoLocalHamiltonianCoefficients);

            oneLocalHamiltonianCoefficients = new Double[] {0, 0 };
            twoLocalHamiltonianCoefficients = new Double[]{ 0, 1, 0, 0};
            var maxCut4 = new ProblemInstance(oneLocalHamiltonianCoefficients, twoLocalHamiltonianCoefficients);

            //END EXAMPLES

            var hybridQaoa = new HybridQaoa(numberOfIterations, p, maxCut4, numberOfRandomStartingPoints, true);

            Solution result = hybridQaoa.RunOptimization();
        }
    }
}