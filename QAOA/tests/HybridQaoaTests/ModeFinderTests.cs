// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.QAOA.HybridQaoaTests
{
    using Microsoft.VisualStudio.TestTools.UnitTesting;
    using System.Collections.Generic;
    using Microsoft.Quantum.QAOA.QaoaHybrid.Helpers;

    [TestClass]
    public class ModeFinderTests
    {
        [TestMethod]
        public void FindModeInBoolListTest()
        {

            var listOfBools = new List<bool[]>
            {
                new[] { false, false, true },
                new[] { false, false, true },
                new[] { false, false, false },
                new[] { false, true, true }
            };

            var expectedResult = new[] {false, false, true};

            var result = ModeFinder.FindModeInBoolList(listOfBools);

            CollectionAssert.AreEqual(expectedResult, result, "Mode bool string not found correctly.");


        }
    }
}

