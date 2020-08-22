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
            bool[] boolsArray1 = { false, false, true };
            bool[] boolsArray2 = { false, false, true };
            bool[] boolsArray3 = { false, false, false };
            bool[] boolsArray4 = { false, true, true };

            var listOfBools = new List<bool[]>
            {
                boolsArray1,
                boolsArray3,
                boolsArray2,
                boolsArray4
            };

            var expectedResult = new[] {false, false, true};

            var result = ModeFinder.FindModeInBoolList(listOfBools);

            CollectionAssert.AreEqual(expectedResult, result, "Mode bool string not found correctly.");


        }
    }
}

