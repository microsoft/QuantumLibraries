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
    //using static FermionTermType.Common;
    //using FermionTerm = FermionTerm;
    //using FermionTermType = FermionTermType;
    using SpinOrbital = SpinOrbital;
    using OrbitalIntegral = OrbitalIntegral;
    using Spin = Spin;
    using Conversion;

    public class OrbitalIntegralToFermionHamiltonianTests
    {

        private static IEnumerable<OrbitalIntegral> orbitalIntegrals = new OrbitalIntegral[]
            {
                    new OrbitalIntegral(new[] { 0,0 }, -1.252477495),
                    new OrbitalIntegral(new[] { 1,1 }, -0.475934275),
                    new OrbitalIntegral(new[] { 0,0,0,0 }, 0.674493166),
                    new OrbitalIntegral(new[] { 0,1,0,1 }, 0.181287518),
                    new OrbitalIntegral(new[] { 0,1,1,0 }, 0.663472101),
                    new OrbitalIntegral(new[] { 1,1,1,1 }, 0.697398010)
            };

        [Fact]
        void BuildOrbitalIntegralHamiltonian()
        {
            var hamiltonian = new OrbitalIntegralHamiltonian();
            hamiltonian.AddTerms(orbitalIntegrals);
        }

        [Fact]
        public void BuildFermionHamiltonian()
        {
            var sourceHamiltonian = new OrbitalIntegralHamiltonian();
            sourceHamiltonian.AddTerms(orbitalIntegrals);

            var targetHamiltonian0 = sourceHamiltonian.ToFermionHamiltonian(SpinOrbital.Config.IndexConvention.Type.HalfUp);

            var targetHamiltonian1 = sourceHamiltonian.ToFermionHamiltonian(SpinOrbital.Config.IndexConvention.Type.UpDown);
        }
        /*
               
            [Fact]
            public void VerifyFermionTermsTest()
            {
                var ham = GenerateTestHamiltonian();
                Assert.True(ham.VerifyFermionTerms());
            }


        }

        public class FermionTermTypeTests
        {

            [Theory]
            [InlineData(true, 1L, new Int64[] { 0, 0, 0, 0 })]
            [InlineData(true, 1L, new Int64[] { 1, 1, 0, 0 })]
            [InlineData(true, 2L, new Int64[] { 1, 0, })]
            [InlineData(true, 2L, new Int64[] { 1, 1, })]
            [InlineData(true, 0L, new Int64[] { })]
            [InlineData(false, 1L, new Int64[] { 0, 0, 1, 0 })]
            [InlineData(false, 1L, new Int64[] { 1, 0, 1, 0 })]
            [InlineData(false, 3L, new Int64[] { 1, 0, })]
            [InlineData(false, 1L, new Int64[] { })]
            public void IsInCanonicalOrderTest(bool pass, Int64 uniqueTerms, Int64[] type)
            {
                var tmp = new FermionTermType(uniqueTerms, type);
                if (pass)
                {
                    Assert.True(tmp.IsInCanonicalOrder());
                }
                else
                {
                    Assert.False(tmp.IsInCanonicalOrder());
                }
            }


            [Fact]
            public void IsInCanonicalOrderCommonTypesTest()
            {
                var tmp = new FermionTermType[] { IdentityTermType,
                PPTermType,
                PQTermType,
                PQQPTermType,
                PQQRTermType,
                PQRSTermType};
                foreach(var item in tmp){
                    Assert.True(item.IsInCanonicalOrder());            
                }
            }
        }
        */
    }
}