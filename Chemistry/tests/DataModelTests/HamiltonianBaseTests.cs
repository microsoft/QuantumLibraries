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
    using static FermionTermType.Common;
    using FermionTerm = FermionTerm;
    using FermionTermType = FermionTermType;
    using SpinOrbital = SpinOrbital;
    using OrbitalIntegral = OrbitalIntegral;
    using Spin = Spin;

    public class FermionHamiltonianTests
    {
       public FermionHamiltonian GenerateTestHamiltonian()
        {
            Int64 nOrbitals = 6;
            Dictionary<FermionTermType, List<FermionTerm>> fermionTerms = new Dictionary<FermionTermType, List<FermionTerm>>
            {
                { PPTermType, new List<FermionTerm>() },
                { PQTermType, new List<FermionTerm>() },
                { PQQPTermType, new List<FermionTerm>() },
                { PQQRTermType, new List<FermionTerm>() },
                { PQRSTermType, new List<FermionTerm>() }
            };

            fermionTerms[PPTermType].Add(new FermionTerm(nOrbitals, new Int64[] { 1, 0 }, new Int64[] { 0, 0 }, 1.0));
            fermionTerms[PPTermType].Add(new FermionTerm(nOrbitals, new Int64[] { 1, 0 }, new Int64[] { 1, 1 }, 1.0));
            fermionTerms[PPTermType].Add(new FermionTerm(nOrbitals, new Int64[] { 1, 0 }, new Int64[] { 2, 2 }, 1.0));

            fermionTerms[PQTermType].Add(new FermionTerm(nOrbitals, new Int64[] { 1, 0 }, new Int64[] { 0, 2 }, 1.0));
            fermionTerms[PQTermType].Add(new FermionTerm(nOrbitals, new Int64[] { 1, 0 }, new Int64[] { 1, 3 }, 1.0));
            fermionTerms[PQTermType].Add(new FermionTerm(nOrbitals, new Int64[] { 1, 0 }, new Int64[] { 2, 6 }, 1.0));

            fermionTerms[PQQPTermType].Add(new FermionTerm(nOrbitals, new Int64[] { 1, 1, 0, 0 }, new Int64[] { 0, 2,2 ,0 }, 1.0));
            fermionTerms[PQQPTermType].Add(new FermionTerm(nOrbitals, new Int64[] { 1, 1, 0, 0 }, new Int64[] { 1, 3, 3,1 }, 1.0));
            fermionTerms[PQQPTermType].Add(new FermionTerm(nOrbitals, new Int64[] { 1, 1, 0, 0 }, new Int64[] { 2, 6, 6, 2 }, 1.0));

            fermionTerms[PQQRTermType].Add(new FermionTerm(nOrbitals, new Int64[] { 1, 1, 0, 0 }, new Int64[] { 0, 2, 2, 1 }, 1.0));
            fermionTerms[PQQRTermType].Add(new FermionTerm(nOrbitals, new Int64[] { 1, 1, 0, 0 }, new Int64[] { 1, 3, 3, 2 }, 1.0));
            fermionTerms[PQQRTermType].Add(new FermionTerm(nOrbitals, new Int64[] { 1, 1, 0, 0 }, new Int64[] { 2, 6, 6, 5 }, 1.0));

            fermionTerms[PQRSTermType].Add(new FermionTerm(nOrbitals, new Int64[] { 1, 1, 0, 0 }, new Int64[] { 0, 2, 4, 3 }, 1.0));
            fermionTerms[PQRSTermType].Add(new FermionTerm(nOrbitals, new Int64[] { 1, 1, 0, 0 }, new Int64[] { 1, 4, 3, 2 }, 1.0));
            fermionTerms[PQRSTermType].Add(new FermionTerm(nOrbitals, new Int64[] { 1, 1, 0, 0 }, new Int64[] { 2, 4,5, 3 }, 1.0));
            return new FermionHamiltonian(fermionTerms, nOrbitals);
        }

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

    public class SpinOrbitalTests
    {
        [Theory]
        [InlineData(10, 1, 0,1)]
        [InlineData(10, 1, 1, 11)]
        [InlineData(10, 0, 1, 10)]
        [InlineData(1, 0, 1, 1)]
        [InlineData(1, 0, 0, 0)]
        public void SpinOrbitalTest(Int64 nOrbitals, Int64 orbital, Int64 spin, Int64 expected){
            Assert.Equal(expected, new SpinOrbital((orbital, spin)).ToInt(nOrbitals));
            Assert.Equal(new SpinOrbital(nOrbitals, expected).ToInt(nOrbitals), expected);
        }
    }

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

   

    public class FermionTermTests
    {

        [Theory]
        [InlineData(3, new Int64[] { 0, 0, 0, 0 }, 1)]
        [InlineData(2, new Int64[] { 0, 1, 3, 0 }, 3)]
        [InlineData(2, new Int64[] { 0, 0, 1, 0 }, 2)]
        [InlineData(3, new Int64[] { 1, 2, 2, 1 }, 2)]
        [InlineData(3, new Int64[] { 4, 3, 2, 1 }, 4)]
        public void UniqueIndicesTests(Int64 norbitals, Int64[] idx, Int64 uniqueIndices)
        {
            var spinOrbitals = idx.Select(o => new SpinOrbital(norbitals, o));
            var coefficient = 1.0;
            var fermionTerm = new FermionTerm(spinOrbitals, coefficient);
            Assert.True(fermionTerm.GetUniqueIndices() == uniqueIndices);
        }

        [Theory]
        [InlineData(true, 10, new Int64[] { 0, 0, 0, 0 }, new Int64[] {0,0,0,0 })]
        [InlineData(true, 10, new Int64[] { 1, 1, 0, 0 }, new Int64[] { 0, 0, 0, 0 })]
        [InlineData(true, 10, new Int64[] { 1, 0, 0, 0 }, new Int64[] { 0, 0, 0, 0 })]
        [InlineData(true, 10, new Int64[] { 1, 1, 0, 0 }, new Int64[] { 1, 2, 2, 1 })]
        [InlineData(true, 10, new Int64[] { 0, 0, 0, 0 }, new Int64[] { 9, 2, 2, 1 })]
        [InlineData(false, 10, new Int64[] { 1, 1, 0, 0 }, new Int64[] { 5, 7, 6, 1 })]
        [InlineData(true, 10, new Int64[] { 1, 1, 1, 0 }, new Int64[] { 1, 19, 19, 3 })]
        [InlineData(true, 10, new Int64[] { 1, 0 }, new Int64[] { 1,5 })]
        [InlineData(true, 10, new Int64[] { 1, 0 }, new Int64[] { 1, 1 })]
        [InlineData(false, 10, new Int64[] { 0, 0, 0, 0 }, new Int64[] { 0, 0, 0, 1 })]
        [InlineData(false, 10, new Int64[] { 1, 1, 0, 0 }, new Int64[] { 0, 1, 1, 2 })]
        [InlineData(false, 10, new Int64[] { 1, 0, 0, 0 }, new Int64[] { 0, 5, 4, 5 })]
        [InlineData(false, 10, new Int64[] { 1, 1, 0, 0 }, new Int64[] { 3, 2, 2, 1 })]
        [InlineData(false, 10, new Int64[] { 0, 0, 0, 0 }, new Int64[] { 1, 2, 2, 3 })]
        [InlineData(false, 10, new Int64[] { 1, 1, 0, 0 }, new Int64[] { 5, 7, 6, 7 })]
        [InlineData(true, 10, new Int64[] { 1, 1, 1, 0 }, new Int64[] { 1, 9, 9, 12 })]
        [InlineData(true, 10, new Int64[] { 1, 1, 0, 0 }, new Int64[] { 0, 1, 2, 0 })]
        [InlineData(false, 10, new Int64[] { 1, 1, 0, 0 }, new Int64[] { 0, 2, 1, 0 })]
        [InlineData(false, 10, new Int64[] { 1, 0 }, new Int64[] { 6, 5 })]
        [InlineData(false, 10, new Int64[] { 1, 1 }, new Int64[] { 6, 5 })]
        public void IsInCanonicalOrderTest(bool pass, Int64 nOrbitals, Int64[] ca, Int64[] idx)
        {
            var tmp = new FermionTerm(nOrbitals, ca, idx, 1.0);
            if (pass)
            {
                Assert.True(tmp.IsInCanonicalOrder());
            }
            else
            {
                Assert.False(tmp.IsInCanonicalOrder());
            }
        }


        [Theory]
        [InlineData(true, 10, new Int64[] { 0, 1 }, new Int64[] { 0, 1 })]
        [InlineData(true, 10, new Int64[] { 0, 1 }, new Int64[] { 1, 0 })]
        [InlineData(true, 10, new Int64[] { 0, 1 }, new Int64[] { 1, 1 })]
        [InlineData(true, 10, new Int64[] { 0, 0, 1, 1 }, new Int64[] { 1, 2, 3, 4 })]
        [InlineData(true, 10, new Int64[] { 1, 1, 0, 0 }, new Int64[] { 0, 1, 2, 0 })]
        [InlineData(true, 10, new Int64[] { 1, 1, 0, 0 }, new Int64[] { 0, 2, 1, 0 })]
        [InlineData(true, 10, new Int64[] { 0, 1, 1, 1 }, new Int64[] { 1, 1, 3, 4 })]
        [InlineData(true, 10, new Int64[] { 0, 1, 1 }, new Int64[] { 0, 0, 1 })]
        [InlineData(true, 10, new Int64[] { 0, 1, 1 }, new Int64[] { 0, 1, 0 })]
        public void ToCanonicalOrderTest(bool pass, Int64 nOrbitals, IEnumerable<Int64> ca, IEnumerable<Int64> idx)
        {
            var tmp = new FermionTerm(nOrbitals, ca.ToArray(), idx.ToArray(), 1.0);
            var newTerms = tmp.ToCanonicalOrder();
            foreach (var newTerm in newTerms)
            {
                Assert.True(newTerm.IsInCanonicalOrder());
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
            foreach (var item in tmp)
            {
                Assert.True(item.IsInCanonicalOrder());
            }
        }

        [Theory]
        [InlineData(true, new Int64[] { 1, 2, 1, 3 }, new Spin[] { Spin.u, Spin.u, Spin.d, Spin.d }, -1.0)]
        [InlineData(true, new Int64[] { 1, 2, 1, 3 }, new Spin[] { Spin.d, Spin.u, Spin.d, Spin.d }, 1.0)]
        public void CreateFermionTermTest(bool pass, Int64[] orbitalIdx, Spin[] spinIdx, Double sign)
        {
            var coeff = 1.0;
            var spinOrbital = orbitalIdx.Zip(spinIdx, (a, b) => new SpinOrbital(a, b));
            var tmp = new FermionTerm(spinOrbital, coeff);
            Assert.True(tmp.IsInCanonicalOrder());
            Assert.True(tmp.coeff == sign);
        }


        [Theory]
        [InlineData(3, new Int64[] {}, 0)]
        [InlineData(3, new Int64[] { 0, 0 }, 1)]
        [InlineData(3, new Int64[] { 0, 1 }, 2)]
        [InlineData(2, new Int64[] { 0, 1, 3, 0 }, 4)]
        [InlineData(3, new Int64[] { 1, 2, 2, 1 }, 3)]
        [InlineData(3, new Int64[] { 4, 3, 2, 1 }, 5)]
        public void GetFermionTermTypeTests(Int64 norbitals, Int64[] idx, Int64 type)
        {
            var tmp = new FermionTermType[] {   IdentityTermType,
                                                PPTermType,
                                                PQTermType,
                                                PQQPTermType,
                                                PQQRTermType,
                                                PQRSTermType};


            var spinOrbitals = idx.Select(o => new SpinOrbital(norbitals, o)).ToArray();
            var coefficient = 1.0;
            var fermionTerm = new FermionTerm(spinOrbitals, coefficient);
            var fermionTermType = fermionTerm.GetFermionTermType();
            Assert.True(fermionTermType == tmp[type]);
        }

    }

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
    }
}