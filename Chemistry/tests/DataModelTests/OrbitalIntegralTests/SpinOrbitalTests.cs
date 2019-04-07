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
    using SpinOrbital = SpinOrbital;
    using Spin = Spin;
  
    public class SpinOrbitalTests
    {
        [Theory]
        [InlineData(10, 1, 0,1)]
        [InlineData(10, 1, 1, 11)]
        [InlineData(10, 0, 1, 10)]
        [InlineData(1, 0, 1, 1)]
        [InlineData(1, 0, 0, 0)]
        public void SpinOrbitalTest(int nOrbitals, int orbital, int spin, int expected){
            Assert.Equal(expected, new SpinOrbital((orbital, spin)).ToInt(SpinOrbital.IndexConvention.HalfUp, nOrbitals));
        }
    }
    
}