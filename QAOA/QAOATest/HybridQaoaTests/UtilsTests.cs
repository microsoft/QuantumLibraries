using Microsoft.VisualStudio.TestTools.UnitTesting;
using System;
using Quantum.QAOA;
using System.Collections.Generic;
using System.Text;

namespace QAOATest.ClassicalOptimizationTests
{
    [TestClass]
    public class UtilsTests
    {
        [TestMethod]
        public void modeOfABoolListTest()
        {
            bool[] boolsArray1 = { false, false, true };
            bool[] boolsArray2 = { false, false, true };
            bool[] boolsArray3 = { false, false, false };
            bool[] boolsArray4 = { false, true, true };

            List<bool[]> listOfBools = new List<bool[]>();
            listOfBools.Add(boolsArray1);
            listOfBools.Add(boolsArray3);
            listOfBools.Add(boolsArray2);
            listOfBools.Add(boolsArray4);

            string expectedResult = "001";

            string result = ClassicalOptimizationUtils.getModeFromBoolList(listOfBools);

            Assert.AreEqual(expectedResult, result, "Mode bool string not found correctly.");


        }

        [TestMethod]
        public void boolStringFromBoolArrayTest()
        {
            bool[] boolsArray = { false, false, true };


            string expectedResult = "001";

            string result = ClassicalOptimizationUtils.getBoolStringFromBoolArray(boolsArray);

            Assert.AreEqual(expectedResult, result, "Bool string not created correctly.");


        }
    }
}

