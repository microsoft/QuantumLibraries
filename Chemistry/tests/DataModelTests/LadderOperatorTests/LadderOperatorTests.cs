// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using Xunit;

using System;
using System.Linq;
using System.Collections.Generic;
using Microsoft.Quantum.Chemistry.LadderOperators;
using Microsoft.Quantum.Chemistry.OrbitalIntegrals;

namespace Microsoft.Quantum.Chemistry.Tests
{

    using static RaisingLowering;
    using LadderSequence = LadderSequence<int>;
    using LadderOperator = LadderOperator<int>;
    using NormalOrderedSequence = NormalOrderedSequence<int>;
    using IndexOrderedSequence = IndexOrderedSequence<int>;

    public class LadderOperatorTests
    {
        [Fact]
        public void Empty()
        {
            var term = new LadderSequence(new List<LadderOperator>());
            var term2 = new LadderSequence(new List<LadderOperator>());

            Dictionary<LadderSequence, double> dictionary = new Dictionary<LadderSequence, double>();
            dictionary.Add(term, 0.5);

            Assert.Equal(0.5, dictionary[term2]);
        }

        [Fact]
        public void NotPassByReference()
        {
            var ints = new[] { 1, 2, 4, 3 };
            var op = new[] { 1, 2, 4, 3 }.ToLadderSequence();
            var newTerm = new LadderSequence(op);
            
            ints[2] = 5;
            Assert.Equal(op.Sequence, newTerm.Sequence);

            op.Sequence[2] = new LadderOperator(RaisingLowering.d, 6);
            Assert.NotEqual(op.Sequence, newTerm.Sequence);

            var newTerm2 = new LadderSequence(ints.ToLadderSequence());
            ints[1] = 0;
            Assert.Equal(2, newTerm2.Sequence[1].Index);

            var op2 = new[] { (u, 2), (d, 5), (d, 8) };
            var newTerm3 = new LadderSequence(op2.ToLadderSequence());

            op2[2] = (d, 4);

            Assert.Equal(8, newTerm3.Sequence[2].Index);

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
            var term = new LadderSequence(idx.ToLadderSequence());
            Assert.True(term.UniqueIndices() == uniqueIndices);
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
            var ladderOperators = ca.Zip(idx, (a, b) => (a == 0 ? RaisingLowering.d : RaisingLowering.u, (int) b)).Select(o => new LadderOperator(o)).ToList();
            var tmp = new LadderSequence(ladderOperators);
            
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
            var ladderOperators = ca.Zip(idx, (a, b) => (a == 0 ? RaisingLowering.d : RaisingLowering.u, (int)b)).ToLadderSequence();
            Assert.Throws<ArgumentException>(() => new IndexOrderedSequence(ladderOperators));
        }
        
    }


    public class NormalOrderedLadderOperatorTests
    {
        [Fact]
        public void NotPassByReference()
        {
            var op = new[] { 1, 2, 4, 3 }.ToLadderSequence();
            var newTerm = new NormalOrderedSequence(op);
            op.Sequence[2] = new LadderOperator(RaisingLowering.d, 5);

            Assert.NotEqual(op.Sequence, newTerm.Sequence);
        }
        
        [Fact]
        public void CommuteToNormalOrder()
        {
            var term = new[] { (d, 2), (u, 5), (d, 8) }.ToLadderSequence();

            var expected = new[] { (u, 5), (d, 2), (d, 8) }.ToLadderSequence();
            var output = term.ToNormalOrder();
            Assert.Equal(expected.Sequence, output.First().Sequence);
            Assert.Equal(-1, output.First().Coefficient);
        }

        [Fact]
        public void ToNormalOrder()
        {
            var term = new[] { (d, 2), (u, 5), (d, 8) }.ToLadderSequence();

            var expected = new[] { (u, 5), (d, 8), (d, 2) }.ToLadderSequence();
            var output = term.ToIndexOrder().ToList();
            Assert.Equal(expected.Sequence, output.First().Sequence);
        }

    }


    public class IndexedlOrderedLadderOperatorTests
    {
        [Fact]
        public void NotPassByReference()
        {
            var op = new[] { 1, 2, 4, 3 }.ToLadderSequence();
            var newTerm = new IndexOrderedSequence(op);
            op.Sequence[2] = new LadderOperator(RaisingLowering.d, 5);

            Assert.NotEqual(op.Sequence, newTerm.Sequence);
        }


        [Fact]
        public void ToIndexOrder()
        {
            var term = new[] { (d, 2), (u, 2), (d, 8) }.ToLadderSequence();

            var expected = new[] {
                new[] { (u, 2), (d, 8), (d, 2) }.ToLadderSequence(1),
                new[] { (d, 8) }.ToLadderSequence(1)
            };
            var output = term.ToIndexOrder();
            Assert.Contains(expected[0].Sequence, output.Select(o => o.Sequence));
            Assert.Contains(expected[1].Sequence, output.Select(o => o.Sequence));
        }


        [Fact]
        public void ToSpinOrbitalIndexOrder()
        {

            var term = new[] { (2, Spin.u), (1, Spin.u) };

            // Test after indexing
            var termInts = term.Select(o => new SpinOrbital(o).ToInt());

            var orderedSpinOrbitals = new IndexOrderedSequence<SpinOrbital>(term.ToLadderSequence());
            var orderedInts = new IndexOrderedSequence<int>(termInts.ToLadderSequence());

        }

    }

}