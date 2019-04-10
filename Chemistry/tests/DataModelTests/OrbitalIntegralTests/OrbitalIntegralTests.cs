// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using Xunit;
using Microsoft.Quantum.Chemistry;
using Microsoft.Quantum.Simulation.Core;

using System.Text.RegularExpressions;
using System;
using System.Linq;
using System.Collections.Generic;

using Microsoft.Quantum.Chemistry.OrbitalIntegrals;

namespace Microsoft.Quantum.Chemistry.Tests
{

    public class OrbitalIntegralTests
    {
        [Theory]
        [InlineData(0, 0, 0, 0, 1)]
        [InlineData(0, 0, 0, 1, 4)]
        [InlineData(0, 0, 1, 0, 4)]
        [InlineData(0, 1, 0, 0, 4)]
        [InlineData(1, 0, 0, 0, 4)]
        [InlineData(1, 1, 0, 0, 4)]
        [InlineData(0, 1, 1, 0, 2)]
        [InlineData(0, 0, 1, 2, 8)]
        [InlineData(0, 1, 2, 0, 4)]
        [InlineData(0, 1, 2, 3, 8)]
        public void OrbitalIntegralEnumerateOrbitalSymmetriesTest(int i, int j, int k, int l, int elements)
        {
            OrbitalIntegral orbitalIntegral = new OrbitalIntegral(new int[] { i, j, k, l });
            var orbitalIntegrals = orbitalIntegral.EnumerateOrbitalSymmetries();
            Assert.Equal(elements, orbitalIntegrals.Length);
        }

        [Theory]
        [InlineData(0, 0, 0, 0, 1)]
        [InlineData(0, 0, 0, 1, 4)]
        [InlineData(0, 0, 1, 0, 4)]
        [InlineData(0, 1, 0, 0, 4)]
        [InlineData(1, 0, 0, 0, 4)]
        [InlineData(1, 1, 0, 0, 4)]
        [InlineData(0, 1, 1, 0, 2)]
        [InlineData(0, 0, 1, 2, 8)]
        [InlineData(0, 1, 2, 0, 4)]
        [InlineData(0, 1, 2, 3, 8)]
        public void OrbitalIntegralEnumerateSpinOrbitalsTest(int i, int j, int k, int l, int elements)
        {
            OrbitalIntegral orbitalIntegral = new OrbitalIntegral(new int[] { i, j, k, l });
            var orbitalIntegrals = orbitalIntegral.EnumerateOrbitalSymmetries();
            var spinOrbitals = orbitalIntegrals.EnumerateSpinOrbitals();
            Assert.Equal(elements * 4, spinOrbitals.Length);
        }

        [Fact]
        public void OrbitalIntegralIndexTest()
        {
            var orb0 = new OrbitalIntegral( 1.0);
            var orb1 = new OrbitalIntegral( 2.0);
            var orb2 = new OrbitalIntegral(new[] { 1, 2, 3, 4 }, 0.5);
            Assert.True(orb0 == orb1);
            Assert.False(orb0 == orb2);
        }
    }
}