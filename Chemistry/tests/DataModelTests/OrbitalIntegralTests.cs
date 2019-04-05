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

    using OrbitalIntegral = OrbitalIntegral;
 
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
        public void OrbitalIntegralEnumerateOrbitalSymmetriesTest(Int64 i, Int64 j, Int64 k, Int64 l, Int64 elements)
        {
            OrbitalIntegral orbitalIntegral = new OrbitalIntegral(new Int64[] { i, j, k, l });
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
        public void OrbitalIntegralEnumerateSpinOrbitalsTest(Int64 i, Int64 j, Int64 k, Int64 l, Int64 elements)
        {
            OrbitalIntegral orbitalIntegral = new OrbitalIntegral(new Int64[] { i, j, k, l });
            var orbitalIntegrals = orbitalIntegral.EnumerateOrbitalSymmetries();
            var spinOrbitals = orbitalIntegrals.EnumerateSpinOrbitals();
            Assert.Equal(elements * 4, spinOrbitals.Length);
        }
    }
}