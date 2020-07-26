using Microsoft.VisualStudio.TestTools.UnitTesting;
using System.Collections.Generic;
using Quantum.QAOA;
using System;
using System.Reflection;

namespace QAOATest.ClassicalOptimizationTests
{

    [TestClass]
    public class ClassicalOptimizationTest
    {

        [TestMethod]
        public void convertDataVectorToVectorsTest()
        {

            ProblemInstance problemInstance = new ProblemInstance(new double[] { 1, 2, 2, -1 }, new double[] { 5, 0, 0, 1, 1, 5, 0, 0, 3, 4, -2, -2, 8, 7, -2, 12 });

            HybridQaoa classicalOptimization = new HybridQaoa(2, 3, problemInstance, 1, new Double[] { 1, 2, 3 }, new Double[] { 4, 5, 6 } );

            ClassicalOptimizationUtils.FreeParamsVector result = ClassicalOptimizationUtils.convertVectorIntoHalves(new Double[] { 1, 2, 3, 4, 5, 6 });

            ClassicalOptimizationUtils.FreeParamsVector dataVectors = new ClassicalOptimizationUtils.FreeParamsVector();
            dataVectors.beta = new double[] { 1, 2, 3 };
            dataVectors.gamma = new double[] { 4, 5, 6 };

            ClassicalOptimizationUtils.FreeParamsVector expectedResult = dataVectors;

            CollectionAssert.AreEqual(expectedResult.beta, result.beta, "Hamiltonian beta value not calculated correctly.");
            CollectionAssert.AreEqual(expectedResult.gamma, result.gamma, "Hamiltonian gamma value not calculated correctly.");

        }

        [TestMethod]
        public void evaluateHamiltonianTest()
        {
            ProblemInstance problemInstance = new ProblemInstance(new double[] { 1, 2, 2, -1 }, new double[] { 5, 0, 0, 1, 1, 5, 0, 0, 3, 4, -2, -2, 8, 7, -2, 12 });

            HybridQaoa classicalOptimization = new HybridQaoa(2, 3, problemInstance, 1, new Double[] { 1, 2, 3 }, new Double[] { 4, 5, 6 });

            Double result = classicalOptimization.evaluateHamiltonian("0011");

            Double expectedResult = -1;

            Assert.AreEqual(expectedResult, result, "Hamiltonian value not calculated correctly.");

        }

        [TestMethod]
        public void evaluateCostFunctionTest()
        {
            ProblemInstance problemInstance = new ProblemInstance(new double[] { 1, 1, 1, 1 }, new double[] { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 });

            HybridQaoa classicalOptimization = new HybridQaoa(2, 1, problemInstance, 1, new double[] { 2 }, new double[] { 3 });


            string optimizationResult = "0101";

            Double result = classicalOptimization.evaluateCostFunction(optimizationResult, new double[] { 5, 3, 2, 1 });

            Double expectedResult = 4;

            Assert.AreEqual(expectedResult, result, "Cost function not calculated correctly.");


        }


    }
}
