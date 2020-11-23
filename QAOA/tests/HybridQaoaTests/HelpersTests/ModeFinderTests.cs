// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

using System;
using System.Linq;

namespace Microsoft.Quantum.Qaoa.HybridQaoaTests
{
    using Microsoft.VisualStudio.TestTools.UnitTesting;
    using System.Collections.Generic;
    using Microsoft.Quantum.Qaoa.QaoaHybrid.Helpers;

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

        [TestMethod]
        public void FindModeInBoolListWithTieTest()
        {

            var listOfBools = new List<bool[]>
            {
                new[] { false, true, true },
                new[] { false, false, true },
                new[] { false, false, false },
                new[] { false, false, true },
                new[] { false, false, false }
            };

            var expectedResult1 = new[] { false, false, true };
            var expectedResult2 = new[] { false, false, false };

            var result = ModeFinder.FindModeInBoolList(listOfBools);

            Assert.IsTrue(expectedResult1.SequenceEqual(result) || expectedResult2.SequenceEqual(result), "Mode bool string not found correctly.");
        }
    }
}

