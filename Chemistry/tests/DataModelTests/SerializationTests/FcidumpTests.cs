// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

#nullable enable

using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

using Microsoft.Quantum.Chemistry.Broombridge;

using Newtonsoft.Json;

using Xunit;

using static Microsoft.Quantum.Chemistry.Tests.Extensions;

namespace Microsoft.Quantum.Chemistry.Tests
{
    public class FciDumpTests : IClassFixture<FciDumpDataFixture>
    {
        FciDumpDataFixture fixture;
        public FciDumpTests(FciDumpDataFixture fixture)
        {
            this.fixture = fixture;
        }

        [Fact]
        public void EnergyOffsetIsZero()
        {
            Assert.Equal(fixture.problem.EnergyOffset.Value, 0.0);
            AssertThat(fixture.problem.EnergyOffset.HasUnits("hartree"));
        }

        [Fact]
        public void CoulombRepulsionIsCorrect()
        {
            Assert.Equal(fixture.problem.CoulombRepulsion.Value, 0.71510433906);
            AssertThat(fixture.problem.CoulombRepulsion.HasUnits("hartree"));
        }

        [Fact]
        public void NumberOfTermsIsCorrect()
        {
            Assert.Equal(fixture.problem.OrbitalIntegralHamiltonian.Terms[TermType.OrbitalIntegral.OneBody].Count, 2);
            Assert.Equal(fixture.problem.OrbitalIntegralHamiltonian.Terms[TermType.OrbitalIntegral.TwoBody].Count, 4);
        }
    }
}
