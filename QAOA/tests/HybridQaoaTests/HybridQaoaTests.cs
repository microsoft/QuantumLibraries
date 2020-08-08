using Microsoft.VisualStudio.TestTools.UnitTesting;
using QAOA.QaoaHybrid;
using System;

namespace Microsoft.Quantum.QAOA.HybridQaoaTests
{

    [TestClass]
    public class ClassicalOptimizationTest
    {

        [TestMethod]
        public void ConvertDataVectorToVectorsTest()
        {

            ProblemInstance problemInstance = new ProblemInstance(new double[] { 1, 2, 2, -1 }, new double[] { 5, 0, 0, 1, 1, 5, 0, 0, 3, 4, -2, -2, 8, 7, -2, 12 });

            HybridQaoa classicalOptimization = new HybridQaoa(2, 3, problemInstance, 1, null, new Double[] { 1, 2, 3 }, new Double[] { 4, 5, 6 } );

            Utils.FreeParamsVector result = Utils.ConvertVectorIntoHalves(new Double[] { 1, 2, 3, 4, 5, 6 });

            Utils.FreeParamsVector dataVectors = new Utils.FreeParamsVector();
            dataVectors.beta = new double[] { 1, 2, 3 };
            dataVectors.gamma = new double[] { 4, 5, 6 };

            Utils.FreeParamsVector expectedResult = dataVectors;

            CollectionAssert.AreEqual(expectedResult.beta, result.beta, "Hamiltonian beta value not calculated correctly.");
            CollectionAssert.AreEqual(expectedResult.gamma, result.gamma, "Hamiltonian gamma value not calculated correctly.");

        }

        [TestMethod]
        public void EvaluateHamiltonianTest()
        {
            ProblemInstance problemInstance = new ProblemInstance(new double[] { 1, 2, 2, -1 }, new double[] { 5, 0, 0, 1, 1, 5, 0, 0, 3, 4, -2, -2, 8, 7, -2, 12 });

            HybridQaoa classicalOptimization = new HybridQaoa(2, 3, problemInstance, 1, null, new Double[] { 1, 2, 3 }, new Double[] { 4, 5, 6 });

            Double result = classicalOptimization.EvaluateHamiltonian("0011");

            Double expectedResult = -1;

            Assert.AreEqual(expectedResult, result, "Hamiltonian value not calculated correctly.");

        }

        [TestMethod]
        public void EvaluateCostFunctionTest()
        {
            ProblemInstance problemInstance = new ProblemInstance(new double[] { 1, 1, 1, 1 }, new double[] { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 });

            HybridQaoa classicalOptimization = new HybridQaoa(2, 1, problemInstance, 1, null, new double[] { 2 }, new double[] { 3 });


            string optimizationResult = "0101";

            Double result = classicalOptimization.EvaluateCostFunction(optimizationResult, new double[] { 5, 3, 2, 1 });

            Double expectedResult = 4;

            Assert.AreEqual(expectedResult, result, "Cost function not calculated correctly.");


        }

        [TestMethod]
        public void RunHybridQaoaTest()
        {
            double[] dh = new Double[] { 0, 0 };
            double[] dJ = new Double[]{ 0, 1,
                               0, 0};

            int numberOfIterations = 50;
            int p = 2;
            int numberOfRandomStartingPoints = 2;

            ProblemInstance simpleMaxCut = new ProblemInstance(dh, dJ);

            HybridQaoa classicalOptimization = new HybridQaoa(numberOfIterations, p, simpleMaxCut, numberOfRandomStartingPoints);
            OptimalSolution optimalSolution = classicalOptimization.RunOptimization();

            string optimizationResult1 = "01";
            string optimizationResult2 = "10";

            string result = optimalSolution.optimalVector;
            Console.WriteLine(result);

            Assert.IsTrue(result.Equals(optimizationResult1) || result.Equals(optimizationResult2), "Hybrid QAOA produced incorrect result.");

        }
    }
}
