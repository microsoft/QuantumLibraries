// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

using Microsoft.VisualStudio.TestTools.UnitTesting;
using System;
using Microsoft.Quantum.QAOA.QaoaHybrid;

namespace Microsoft.Quantum.QAOA.HybridQaoaTests
{

    [TestClass]
    public class ClassicalOptimizationTest
    {

        [TestMethod]
        public void ConvertDataVectorToVectorsTest()
        {

            var result = Utils.ConvertVectorIntoHalves(new Double[] { 1, 2, 3, 4, 5, 6 });

            var dataVectors = new Utils.FreeParameters
            {
                Beta = new double[] { 1, 2, 3 },
                Gamma = new double[] { 4, 5, 6 }
            };

            var expectedResult = dataVectors;

            CollectionAssert.AreEqual(expectedResult.Beta, result.Beta, "Hamiltonian beta value not calculated correctly.");
            CollectionAssert.AreEqual(expectedResult.Gamma, result.Gamma, "Hamiltonian gamma value not calculated correctly.");

        }

        
        [TestMethod]
        public void EvaluateCostFunctionTest()
        {
            ProblemInstance problemInstance = new ProblemInstance(new double[] { 1, 1, 1, 1 }, new double[] { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 });

            HybridQaoa classicalOptimization = new HybridQaoa(2, 1, problemInstance, 1, false, new double[] { 2 }, new double[] { 3 });

            var optimizationResult = new []{false, true, false, true};

            var result = classicalOptimization.EvaluateCostFunction(optimizationResult, new double[] { 5, 3, 2, 1 });

            var expectedResult = 4;

            Assert.AreEqual(expectedResult, result, "Cost function not calculated correctly.");


        }

        [TestMethod]
        public void RunHybridQaoaTest()
        {
            var dh = new double[] { 0, 0 };
            var dJ = new double[]{ 0, 1, 0, 0};

            var numberOfIterations = 50;
            var p = 2;
            var numberOfRandomStartingPoints = 2;

            var simpleMaxCut = new ProblemInstance(dh, dJ);

            var classicalOptimization = new HybridQaoa(numberOfIterations, p, simpleMaxCut, numberOfRandomStartingPoints);
            var optimalSolution = classicalOptimization.RunOptimization();

            var optimizationResult1 = new[] {false, true};
            var optimizationResult2 = new[] {true, false};

            var result = optimalSolution.SolutionVector;

            if (result[0] == false)
            {
                CollectionAssert.AreEqual(result, optimizationResult1, "Hybrid QAOA produced incorrect result.");
            }
            else
            {
                CollectionAssert.AreEqual(result, optimizationResult2, "Hybrid QAOA produced incorrect result.");
            }
        }
    }
}
