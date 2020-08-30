// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Qaoa.HybridQaoaTests
{
    using Microsoft.VisualStudio.TestTools.UnitTesting;
    using Microsoft.Quantum.Qaoa.QaoaHybrid;

    [TestClass]
    public class ClassicalOptimizationTest
    {

        [TestMethod]
        public void EvaluateCostFunctionTest()
        {
            QaoaProblemInstance qaoaProblemInstance = new QaoaProblemInstance(new double[] { 1, 1, 1, 1 }, new double[] { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 });

            HybridQaoa classicalOptimization = new HybridQaoa(2, 1, qaoaProblemInstance);

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
            var NHamiltonianApplications = 2;
            var numberOfRandomStartingPoints = 2;

            var simpleMaxCut = new QaoaProblemInstance(dh, dJ);

            var classicalOptimization = new HybridQaoa(numberOfIterations, NHamiltonianApplications, simpleMaxCut);
            var optimalSolution = classicalOptimization.RunOptimization(numberOfRandomStartingPoints);

            var optimizationResult1 = new[] {false, true};
            var optimizationResult2 = new[] {true, false};

            var result = optimalSolution.SolutionVector;

            CollectionAssert.AreEqual(result, result[0] ? optimizationResult2 : optimizationResult1, "Hybrid QAOA with random parameters produced incorrect result.");
        }

        [TestMethod]
        public void RunOptimizationUserParamsTest()
        {
            var dh = new double[] { 0, 0 };
            var dJ = new double[] { 0, 1, 0, 0 };

            var numberOfIterations = 60;
            var NHamiltonianApplications = 3;

            var simpleMaxCut = new QaoaProblemInstance(dh, dJ);

            var initialBeta = new double[] { 0, 0, 0 };
            var initialGamma = new double[] { 0, 0, 0 };

            var qaoaParameters = new QaoaParameters(initialBeta, initialGamma);

            var classicalOptimization = new HybridQaoa(numberOfIterations, NHamiltonianApplications, simpleMaxCut);
            var optimalSolution = classicalOptimization.RunOptimization(qaoaParameters);

            var optimizationResult1 = new[] { false, true };
            var optimizationResult2 = new[] { true, false };

            var result = optimalSolution.SolutionVector;

            CollectionAssert.AreEqual(result, result[0] ? optimizationResult2 : optimizationResult1, "Hybrid QAOA with random parameters produced incorrect result.");
        }
    }
}
