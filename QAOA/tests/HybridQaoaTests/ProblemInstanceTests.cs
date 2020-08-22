// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.QAOA.HybridQaoaTests
{
    using Microsoft.Quantum.QAOA.QaoaHybrid;
    using Microsoft.VisualStudio.TestTools.UnitTesting;

    [TestClass]
    class ProblemInstanceTests
    {

        
        [TestMethod]
        public void EvaluateHamiltonianTest()
        {
            var problemInstance = new ProblemInstance(new double[] { 1, 2, 2, -1 }, new double[] { 5, 0, 0, 1, 1, 5, 0, 0, 3, 4, -2, -2, 8, 7, -2, 12 });

           var result = problemInstance.EvaluateHamiltonian(new [] {false,false,true,true});

           const int expectedResult = -1;

           Assert.AreEqual(expectedResult, result, "Hamiltonian value not calculated correctly.");


        }
        
}
}
