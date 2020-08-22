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
        public void EvaluateCostFunctionTest()
        {
            ProblemInstance problemInstance = new ProblemInstance(new double[] { 1, 1, 1, 1 }, new double[] { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 });

            HybridQaoa classicalOptimization = new HybridQaoa(2, 1, problemInstance);

            var optimizationResult = new []{false, true, false, true};

            var result = classicalOptimization.EvaluateCostFunction(optimizationResult, new double[] { 5, 3, 2, 1 });

            var expectedResult = 4;

            Assert.AreEqual(expectedResult, result, "Cost function not calculated correctly.");


        }

        [TestMethod]
        public void RunOptimizationRandomParamsTest()
        {
            var dh = new double[] { 0, 0 };
            var dJ = new double[]{ 0, 1, 0, 0};

            var numberOfIterations = 50;
            var p = 2;
            var numberOfRandomStartingPoints = 2;

            var simpleMaxCut = new ProblemInstance(dh, dJ);

            var classicalOptimization = new HybridQaoa(numberOfIterations, p, simpleMaxCut);
            var optimalSolution = classicalOptimization.RunOptimization(numberOfRandomStartingPoints);

            var optimizationResult1 = new[] {false, true};
            var optimizationResult2 = new[] {true, false};

            var result = optimalSolution.SolutionVector;

            if (result[0] == false)
            {
                CollectionAssert.AreEqual(result, optimizationResult1, "Hybrid QAOA with random parameters produced incorrect result.");
            }
            else
            {
                CollectionAssert.AreEqual(result, optimizationResult2, "Hybrid QAOA with random parameters produced incorrect result.");
            }
        }

        [TestMethod]
        public void RunOptimizationUserParamsTest()
        {
            var dh = new double[] { 0, 0 };
            var dJ = new double[] { 0, 1, 0, 0 };

            var numberOfIterations = 60;
            var p = 3;

            var simpleMaxCut = new ProblemInstance(dh, dJ);

            var initialBeta = new double[] { 0, 0, 0 };
            var initialGamma = new double[] { 0, 0, 0 };

            var qaoaParameters = new QaoaParameters(initialBeta, initialGamma);

            var classicalOptimization = new HybridQaoa(numberOfIterations, p, simpleMaxCut);
            var optimalSolution = classicalOptimization.RunOptimization(qaoaParameters);

            var optimizationResult1 = new[] { false, true };
            var optimizationResult2 = new[] { true, false };

            var result = optimalSolution.SolutionVector;

            if (result[0] == false)
            {
                CollectionAssert.AreEqual(result, optimizationResult1, "Hybrid QAOA with user's parameters produced incorrect result.");
            }
            else
            {
                CollectionAssert.AreEqual(result, optimizationResult2, "Hybrid QAOA with user's parameters produced incorrect result.");
            }
        }
    }
}
