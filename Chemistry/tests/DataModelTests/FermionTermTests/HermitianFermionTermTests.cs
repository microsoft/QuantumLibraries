// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using Xunit;
using Microsoft.Quantum.Chemistry;
using Microsoft.Quantum.Simulation.Core;

using System.Text.RegularExpressions;
using System;
using System.Linq;
using System.Collections.Generic;

using Microsoft.Quantum.Chemistry.LadderOperators;
using Microsoft.Quantum.Chemistry.Fermion;
using Microsoft.Quantum.Chemistry.OrbitalIntegrals;

namespace Microsoft.Quantum.Chemistry.Tests
{
    using FermionTermHermitian = FermionTermHermitian;
    using static TermType.Fermion;
    using Op = LadderSequence;
    using static LadderType;

    public class HermitianFermionTermTests
    {
        [Fact]
        public void EmptyHermitianFermionTerm()
        {
            var term = new FermionTermHermitian(new LadderSequence());
            var term2 = new FermionTermHermitian(new LadderSequence());

            Dictionary<FermionTermHermitian, double> dictionary = new Dictionary<FermionTermHermitian, double>();
            dictionary.Add(term, 0.5);

            Assert.Equal(0.5, dictionary[term2]);
        }

        [Fact]
        public void Inheritance()
        {
            var term = new NormalOrderedLadderSequence(new LadderSequence());
            var term2 = new FermionTermHermitian(term);
        }


        [Theory]
        [InlineData(3, new[] { 0, 0, 0, 0 }, 1)]
        [InlineData(2, new[] { 0, 1, 3, 0 }, 3)]
        [InlineData(2, new[] { 0, 0, 1, 0 }, 2)]
        [InlineData(3, new[] { 1, 2, 2, 1 }, 2)]
        [InlineData(3, new[] { 4, 3, 2, 1 }, 4)]
        public void UniqueIndicesTests(int norbitals, int[] idx, int uniqueIndices)
        {
            var coefficient = 1.0;
            var HermitianFermionTerm = new FermionTermHermitian(idx.ToLadderSequence());
            Assert.True(HermitianFermionTerm.GetUniqueIndices() == uniqueIndices);
        }

        [Theory]
        [InlineData(true, 10, new int[] { 0, 0, 0, 0 }, new int[] {0,0,0,0 })]
        [InlineData(true, 10, new int[] { 1, 1, 0, 0 }, new int[] { 0, 0, 0, 0 })]
        [InlineData(true, 10, new int[] { 1, 0, 0, 0 }, new int[] { 0, 0, 0, 0 })]
        [InlineData(true, 10, new int[] { 1, 1, 0, 0 }, new int[] { 1, 2, 2, 1 })]
        [InlineData(true, 10, new int[] { 0, 0, 0, 0 }, new int[] { 9, 2, 2, 1 })]
        [InlineData(false, 10, new int[] { 1, 1, 0, 0 }, new int[] { 5, 7, 6, 1 })]
        [InlineData(true, 10, new int[] { 1, 1, 1, 0 }, new int[] { 1, 19, 19, 3 })]
        [InlineData(true, 10, new int[] { 1, 0 }, new int[] { 1,5 })]
        [InlineData(true, 10, new int[] { 1, 0 }, new int[] { 1, 1 })]
        [InlineData(false, 10, new int[] { 0, 0, 0, 0 }, new int[] { 0, 0, 0, 1 })]
        [InlineData(false, 10, new int[] { 1, 1, 0, 0 }, new int[] { 0, 1, 1, 2 })]
        [InlineData(false, 10, new int[] { 1, 0, 0, 0 }, new int[] { 0, 5, 4, 5 })]
        [InlineData(false, 10, new int[] { 1, 1, 0, 0 }, new int[] { 3, 2, 2, 1 })]
        [InlineData(false, 10, new int[] { 0, 0, 0, 0 }, new int[] { 1, 2, 2, 3 })]
        [InlineData(false, 10, new int[] { 1, 1, 0, 0 }, new int[] { 5, 7, 6, 7 })]
        [InlineData(true, 10, new int[] { 1, 1, 1, 0 }, new int[] { 1, 9, 9, 12 })]
        [InlineData(true, 10, new int[] { 1, 1, 0, 0 }, new int[] { 0, 1, 2, 0 })]
        [InlineData(false, 10, new int[] { 1, 1, 0, 0 }, new int[] { 0, 2, 1, 0 })]
        [InlineData(false, 10, new int[] { 1, 0 }, new int[] { 6, 5 })]
        [InlineData(false, 10, new int[] { 1, 1 }, new int[] { 6, 5 })]
        public void IsInCanonicalOrderTest(bool pass, int nOrbitals, int[] ca, int[] idx)
        {
            var ladderOperators = ca.Zip(idx, (a, b) => (a == 0 ? LadderType.d : LadderType.u, (int)b)).ToLadderSequence();
            var tmp = new FermionTermHermitian(ladderOperators);

            Assert.True(tmp.IsInCanonicalOrder());
        }

        [Theory]
        [InlineData(true, 10, new int[] { 0, 1 }, new int[] { 0, 1 })]
        [InlineData(true, 10, new int[] { 0, 1 }, new int[] { 1, 0 })]
        [InlineData(true, 10, new int[] { 0, 0, 1, 1 }, new int[] { 1, 2, 3, 4 })]
        [InlineData(true, 10, new int[] { 0, 1 }, new int[] { 1, 1 })]
        [InlineData(true, 10, new int[] { 0, 1, 1, 1 }, new int[] { 1, 1, 3, 4 })]
        [InlineData(true, 10, new int[] { 0, 1, 1 }, new int[] { 0, 0, 1 })]
        [InlineData(true, 10, new int[] { 0, 1, 1 }, new int[] { 0, 1, 0 })]
        public void NotNormalOrderedTest(bool pass, int nOrbitals, IEnumerable<int> ca, IEnumerable<int> idx)
        {
            var ladderOperators = ca.Zip(idx, (a, b) => (a == 0 ? LadderType.d : LadderType.u, (int)b)).ToLadderSequence();
            Assert.Throws<ArgumentException>(() => new FermionTermHermitian(ladderOperators));
        }

        
        /*
        [Theory]
        [InlineData(true, new int[] { 1, 2, 1, 3 }, new Spin[] { Spin.u, Spin.u, Spin.d, Spin.d }, -1.0)]
        [InlineData(true, new int[] { 1, 2, 1, 3 }, new Spin[] { Spin.d, Spin.u, Spin.d, Spin.d }, 1.0)]
        public void CreateHermitianFermionTermTest(bool pass, int[] orbitalIdx, Spin[] spinIdx, Double sign)
        {
            var coeff = 1.0;
            var spinOrbital = orbitalIdx.Zip(spinIdx, (a, b) => new SpinOrbital(a, b));
            var tmp = new FermionTermHermitian((IEnumerable<SpinOrbital>)spinOrbital, (double)coeff);
            Assert.True(tmp.IsInCanonicalOrder());
            Assert.True(tmp.coeff == sign);
        }
        */
        
        [Theory]
        [InlineData(3, new int[] {}, 0)]
        [InlineData(3, new int[] { 0, 0 }, 1)]
        [InlineData(3, new int[] { 0, 1 }, 2)]
        [InlineData(2, new int[] { 0, 1, 3, 0 }, 4)]
        [InlineData(3, new int[] { 1, 2, 2, 1 }, 3)]
        [InlineData(3, new int[] { 4, 3, 2, 1 }, 5)]
        public void GetFermionTermTypeTests(int norbitals, int[] idx, int type)
        {
            var tmp = new TermType.Fermion[] {  Identity,
                                                PP,
                                                PQ,
                                                PQQP,
                                                PQQR,
                                                PQRS};


            var fermionTerm = new FermionTermHermitian(idx.ToLadderSequence());
            var fermionTermType = fermionTerm.GetTermType();
            Assert.True(fermionTermType == tmp[type]);
        }

        [Fact]
        public void EqualityTest()
        {
            var term0 = new FermionTermHermitian(new[] { (u, 0), (u, 0) }.ToLadderSequence());
            var term1 = new FermionTermHermitian(new[] { (u, 0), (u, 0) }.ToLadderSequence());

            Assert.True(term0 == term1);
            Assert.Equal(term0, term1);
        }

        public (LadderSequence, LadderSequence) CanonicalOrderCorrectnessHelper(int test)
        {
            switch (test)
            {
                case 0:
                    return (new[] { (u, 0), (u, 0) }.ToLadderSequence(),
                            new[] { (u, 0), (u, 0) }.ToLadderSequence());
                case 1:
                    return (new[] { (u, 0), (u, 1) }.ToLadderSequence(-1),
                            new[] { (u, 1), (u, 0) }.ToLadderSequence());
                case 2:
                    return (new[] { (u, 0), (u, 1) }.ToLadderSequence(-1),
                            new[] { (d, 0), (d, 1) }.ToLadderSequence());
                case 3:
                    return (new[] { (u, 0), (d, 1) }.ToLadderSequence(),
                            new[] { (u, 1), (d, 0) }.ToLadderSequence());
                case 4:
                    return (new[] { (u, 2), (u, 5), (u, 9), (d, 10), (d, 5) }.ToLadderSequence(),
                            new[] { (u, 5), (u, 2), (u, 9), (d, 5), (d, 10) }.ToLadderSequence());
                case 5:
                    return (new[] { (u, 5), (u, 9), (u, 10), (d, 5), (d, 2) }.ToLadderSequence(-1),
                            new[] { (u, 5), (u, 2), (d, 9), (d, 5), (d, 10) }.ToLadderSequence());
                default:
                    return (new int[] { }.ToLadderSequence(), new int[] { }.ToLadderSequence());
            }
        }


        [Theory]
        [InlineData(0)]
        [InlineData(1)]
        [InlineData(2)]
        [InlineData(3)]
        [InlineData(4)]
        [InlineData(5)]
        [InlineData(6)]
        public void CanonicalCorrectness(int test)
        {
            var (expected, input) = CanonicalOrderCorrectnessHelper(test);
            var term = new FermionTermHermitian(input);

            Assert.Equal(expected.sequence, term.sequence);
            Assert.Equal(expected.coefficient, term.coefficient);
        }

        /*
[Theory]
[InlineData(true, new int[] { 1, 2, 1, 3 }, new Spin[] { Spin.u, Spin.u, Spin.d, Spin.d }, -1.0)]
[InlineData(true, new int[] { 1, 2, 1, 3 }, new Spin[] { Spin.d, Spin.u, Spin.d, Spin.d }, 1.0)]
public void CreateFermionTermTest(bool pass, int[] orbitalIdx, Spin[] spinIdx, Double sign)
{
    var coeff = 1.0;
    var spinOrbital = orbitalIdx.Zip(spinIdx, (a, b) => new SpinOrbital(a, b));
    var tmp = new FermionTerm((IEnumerable<SpinOrbital>)spinOrbital, (double)coeff);
    Assert.True(tmp.IsInCanonicalOrder());
    Assert.True(tmp.coeff == sign);
}
*/
    }

}