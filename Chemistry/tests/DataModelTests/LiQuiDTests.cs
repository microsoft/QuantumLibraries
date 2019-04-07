// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using Xunit;
using Microsoft.Quantum.Chemistry;
using Microsoft.Quantum.Simulation.Core;
using System.Text.RegularExpressions;
using System;
using System.Linq;
using System.Collections.Generic;
using System.IO;

using YamlDotNet.Core;
using YamlDotNet.Serialization;

using Microsoft.Quantum.Chemistry.OrbitalIntegrals;

namespace Microsoft.Quantum.Chemistry.Tests
{
    using static TermType.OrbitalIntegral;
    using static OrbitalIntegral.Convention;


    public class LoadFromLiquidTests
    {
        [Theory]
        [MemberData(nameof(LiquidOrbitalsData))]
        public void LoadFromLiquidTest(string line, TermType.OrbitalIntegral termType, OrbitalIntegral term)
        {
            //var test = terms.Item1;
            //string[] lines, FermionTermType termType, FermionTerm[] terms
            var hamiltonian = LiQuiD.LoadFromLiquid(line);
            Assert.True(hamiltonian.terms.ContainsKey(termType));
            // Check that expected terms are found

            var orb = hamiltonian.terms[termType].Keys.First();
            // Check index
            //Assert.True(term == orb);

            // Check coefficient
            Assert.Equal(term.Coefficient, hamiltonian.terms[termType][term.ToCanonicalForm()].Value);
        }

        public static IEnumerable<object[]> LiquidOrbitalsData =>
            new List<object[]>
            {
                new object[] { "0,0=1.0", OneBody, new OrbitalIntegral(new[] {0,0},1.0, Dirac) },
                new object[] { "0,1=1.0", OneBody, new OrbitalIntegral(new[] {0,1},1.0, Dirac) },
                new object[] { "0,1,1,0=1.0", TwoBody, new OrbitalIntegral(new[] {0,1,1,0},1.0, Dirac) },
                new object[] { "0,1,0,1=1.0", TwoBody, new OrbitalIntegral(new[] {0,1,0,1},1.0, Dirac) },
                new object[] { "0,1,0,1=1.0", TwoBody, new OrbitalIntegral(new[] {0,1,0,1},1.0, Dirac) },
                new object[] { "0,1,0,0=1.0", TwoBody, new OrbitalIntegral(new[] {0,1,0,0},1.0, Dirac) },
                new object[] { "0,0,1,2=1.0", TwoBody, new OrbitalIntegral(new[] {0,0,1,2},1.0, Dirac) },
                new object[] { "0,0,1,2=1.0", TwoBody, new OrbitalIntegral(new[] {0,0,1,2},1.0, Dirac) },
                new object[] { "0,1,2,0=1.0", TwoBody, new OrbitalIntegral(new[] {0,1,2,0},1.0, Dirac) },
                new object[] { "0,1,2,3=1.0", TwoBody, new OrbitalIntegral(new[] {0,1,2,3},1.0, Dirac) },
                new object[] { "3,1,0,2=1.0", TwoBody, new OrbitalIntegral(new[] {3,1,0,2},1.0, Dirac) }
            };


    }


}