// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Qaoa.HybridQaoaTests
{
    using Microsoft.VisualStudio.TestTools.UnitTesting;
    using Qaoa.QaoaHybrid.Helpers;

    [TestClass]
    class ArrayToStringConverterTests
    {
        [TestMethod]
        public void ConvertBoolArrayToStringTest()
        {
            var boolsArray = new[] { false, false, true };

            var expectedResult = "001";

            var result = ArrayToStringConverter.ConvertBoolArrayToString(boolsArray);

            Assert.AreEqual(expectedResult, result, "Bool string not created correctly.");


        }

    }
}
