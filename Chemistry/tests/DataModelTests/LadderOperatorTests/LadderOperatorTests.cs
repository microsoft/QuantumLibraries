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
    using HermitianFermionTerm = HermitianFermionTerm;
    using SpinOrbital = SpinOrbital;
    using static TermType.Fermion;

    using static LadderOperator.Type;

    public class LadderOperatorTests
    {
        [Fact]
        public void Empty()
        {
            var term = new LadderOperatorSequence(new List<LadderOperator>());
            var term2 = new LadderOperatorSequence(new List<LadderOperator>());

            Dictionary<LadderOperatorSequence, double> dictionary = new Dictionary<LadderOperatorSequence, double>();
            dictionary.Add(term, 0.5);

            Assert.Equal(0.5, dictionary[term2]);
        }

        [Fact]
        public void NotPassByReference()
        {
            var ints = new[] { 1, 2, 4, 3 };
            var op = new LadderOperatorSequence(new[] { 1, 2, 4, 3 });
            var newTerm = new LadderOperatorSequence(op);
            
            ints[2] = 5;
            Assert.Equal(op.sequence, newTerm.sequence);

            op.sequence[2] = new LadderOperator(LadderOperator.Type.d, 6);
            Assert.NotEqual(op.sequence, newTerm.sequence);

            var newTerm2 = new LadderOperatorSequence(ints);
            ints[1] = 0;
            Assert.Equal(2, newTerm2.sequence[1].index);

            var op2 = new[] { (u, 2), (d, 5), (d, 8) };
            var newTerm3 = new LadderOperatorSequence(op2);

            op2[2] = (d, 4);

            Assert.Equal(8, newTerm3.sequence[2].index);

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
            var term = new LadderOperatorSequence(idx);
            Assert.True(term.GetUniqueIndices() == uniqueIndices);
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
            var ladderOperators = ca.Zip(idx, (a, b) => (a == 0 ? LadderOperator.Type.d : LadderOperator.Type.u, (int) b)).Select(o => new LadderOperator(o)).ToList();
            var tmp = new LadderOperatorSequence(ladderOperators);
            
                Assert.True(tmp.IsInNormalOrder());
            
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
            var ladderOperators = ca.Zip(idx, (a, b) => (a == 0 ? LadderOperator.Type.d : LadderOperator.Type.u, (int)b)).Select(o => new LadderOperator(o)).ToList();
            Assert.Throws<ArgumentException>(() => new IndexOrderedLadderOperators(ladderOperators));
        }
        
    }


    public class NormalOrderedLadderOperatorTests
    {
        [Fact]
        public void NotPassByReference()
        {
            var op = new LadderOperatorSequence(new[] { 1, 2, 4, 3 });
            var newTerm = new NormalOrderedLadderOperators(op);
            op.sequence[2] = new LadderOperator(LadderOperator.Type.d, 5);

            Assert.NotEqual(op.sequence, newTerm.sequence);
        }
        
        [Fact]
        public void CommuteToNormalOrder()
        {
            var term = new LadderOperatorSequence(new[] { (d, 2), (u, 5), (d, 8) });

            var expected = new LadderOperatorSequence(new[] { (u, 5), (d, 2), (d, 8) });
            var output = NormalOrderedLadderOperators.CreateNormalOrder(term).ToList();
            Assert.Equal(expected.sequence, output.First().sequence);
            Assert.Equal(-1, output.First().coefficient);
        }

        [Fact]
        public void CreateNormalOrder()
        {
            var term = new LadderOperatorSequence(new[] { (d, 2), (u, 5), (d, 8) });

            var expected = new LadderOperatorSequence(new[] { (u, 5), (d, 8), (d, 2) });
            var output = IndexOrderedLadderOperators.CreateIndexOrder(term).ToList();
            Assert.Equal(expected.sequence, output.First().sequence);
        }

    }


    public class IndexedlOrderedLadderOperatorTests
    {
        [Fact]
        public void NotPassByReference()
        {
            var op = new LadderOperatorSequence(new[] { 1, 2, 4, 3 });
            var newTerm = new IndexOrderedLadderOperators(op);
            op.sequence[2] = new LadderOperator(LadderOperator.Type.d, 5);

            Assert.NotEqual(op.sequence, newTerm.sequence);
        }


        [Fact]
        public void CreateIndexOrder()
        {
            var term = new LadderOperatorSequence(new[] { (d, 2), (u, 2), (d, 8) });

            var expected = new[] {
                new LadderOperatorSequence(new[] { (u, 2), (d, 8), (d, 2) }, 1),
                new LadderOperatorSequence(new[] { (d, 8) })
            };
            var output = IndexOrderedLadderOperators.CreateIndexOrder(term).ToList();
            Assert.Contains(expected[0].sequence, output.Select(o => o.sequence));
            Assert.Contains(expected[1].sequence, output.Select(o => o.sequence));
        }

    }

}