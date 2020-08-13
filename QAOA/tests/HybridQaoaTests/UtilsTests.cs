// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

using Microsoft.VisualStudio.TestTools.UnitTesting;
using System.Collections.Generic;
using Microsoft.Quantum.QAOA.QaoaHybrid;

namespace Microsoft.Quantum.QAOA.HybridQaoaTests
{
    [TestClass]
    public class UtilsTests
    {
        [TestMethod]
        public void ModeOfABoolListTest()
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

            var result = Utils.GetModeFromBoolList(listOfBools);

            CollectionAssert.AreEqual(expectedResult, result, "Mode bool string not found correctly.");


        }

        [TestMethod]
        public void BoolStringFromBoolArrayTest()
        {
            var boolsArray = new [] { false, false, true };

            var expectedResult = "001";

            var result = Utils.GetBoolStringFromBoolArray(boolsArray);

            Assert.AreEqual(expectedResult, result, "Bool string not created correctly.");


        }
    }
}

