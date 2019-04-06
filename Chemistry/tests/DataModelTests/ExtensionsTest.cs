// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using Xunit;
using Microsoft.Quantum.Chemistry;
using Microsoft.Quantum.Simulation.Core;

using System.Text.RegularExpressions;
using System;
using System.Linq;
using System.Collections.Generic;

namespace Microsoft.Quantum.Chemistry.Tests
{
   public class ComparerTests
    {
        [Fact]
        public void IsIntArrayAscendingTest()
        {
            Assert.True((new Int64[] {  }).IsIntArrayAscending());
            Assert.True((new Int64[] { 1 }).IsIntArrayAscending());
            Assert.True((new Int64[] { 1, 2, 3, 4 }).IsIntArrayAscending());
            Assert.True((new Int64[] { 1, 1, 3, 3 }).IsIntArrayAscending());
            Assert.False((new Int64[] { 1, 1, 0, 1 }).IsIntArrayAscending());
            Assert.False((new Int64[] { 10, 9, 7, 4 }).IsIntArrayAscending());
            Assert.False((new Int64[] { 1, 2, 3, 2 }).IsIntArrayAscending());
        }

        [Fact]
        public void CompareIntArrayTest()
        {
            Assert.Equal(0, Extensions.CompareIntArray(new Int64[] { 1, 2, 3, 4 }, new Int64[] { 1, 2, 3, 4 }));
            Assert.Equal(1, Extensions.CompareIntArray(new Int64[] { 1, 2, 3, 4 }, new Int64[] { 1, 2, 3, 3 }));
            Assert.Equal(-1, Extensions.CompareIntArray(new Int64[] { 1, 2, 3, 3 }, new Int64[] { 1, 2, 3, 4 }));
            Assert.Equal(0, Extensions.CompareIntArray(new Int64[] { 2, 2, 2, 2 }, new Int64[] { 2, 2, 2, 2 }));
            Assert.Equal(1, Extensions.CompareIntArray(new Int64[] { 2, 2, 3, 4 }, new Int64[] { 1, 2, 10, 3 }));
            Assert.Equal(-1, Extensions.CompareIntArray(new Int64[] { 1, 3, 3, 3 }, new Int64[] { 10, 2, 3, 4 }));
        }

        [Fact]
        public void PrintIntArray()
        {
            var arr = new[] { 1, 2, 3, 4 };
            var arrString = arr.ToString();

            var arrHashSet = new HashSet<int>(arr) ;
            var arrHashSetString = arr.ToString();

        }
    }
}