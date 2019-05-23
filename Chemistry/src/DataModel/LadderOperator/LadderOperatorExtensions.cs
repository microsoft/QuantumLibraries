// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using System;
using System.Linq;
using System.Collections.Generic;

namespace Microsoft.Quantum.Chemistry.LadderOperators
{ 
    public static partial class Extensions
    {
        #region Convenience constructors

        /// <summary>
        /// Construct a sequence of ladder operators from sequence of tuples each
        /// specifying whether it is a raising or lowering term, and its index.
        /// </summary>
        /// <param name="setSequence">Sequence of ladder operators.</param>
        /// <param name="setSign">Set the sign coefficient of the sequence.</param>
        /// <returns>
        /// Sequence of ladder operators.
        /// </returns>
        /// <example>
        /// // Construct a sequence a ladder operators 1^ 2^ 3 4
        /// var tmp = new[] { (u, 1), (u, 2), (d, 3), (d, 4) }.ToLadderSequence();
        /// </example>
        public static LadderSequence<int> ToLadderSequence(this IEnumerable<(RaisingLowering, int)> setSequence, int setSign = 1) =>
            new LadderSequence<int>(setSequence.Select(o => new LadderOperator<int>(o)), setSign);

        /// <summary>
        /// Construct a sequence of ladder operators from an even-length sequence of integers.
        /// </summary>
        /// <param name="indices">Even-length sequence of integers.</param>
        /// <returns>
        /// Sequence of ladder operators with an equal number of creation and annihilation terms
        /// that are normal-ordered.
        /// </returns>
        /// <example>
        /// <code>
        /// // The following two return the same ladder operator sequence.
        /// var seq = new[] { 1, 2, 3, 4 }.ToLadderSequence();
        /// var expected = new[] { (u, 1), (u, 2), (d, 3), (d, 4) }.ToLadderSequence();
        /// </code>
        /// </example>
        public static LadderSequence<int> ToLadderSequence(this IEnumerable<int> indices) => indices.ToArray();
        #endregion

        #region Reordering methods

        /// <summary>
        ///  Converts a <see cref="LadderSequence"/> to normal order. 
        ///  In general, this can generate new terms and modifies the coefficient.
        /// </summary>
        public static HashSet<NormalOrderedSequence<TIndex>> ToNormalOrder<TIndex>(this LadderSequence<TIndex> ladderOperator)
            where TIndex : IEquatable<TIndex>
        {
            // Recursively anti-commute creation to the left.
            var TmpTerms = new Stack<LadderSequence<TIndex>>();
            var NewTerms = new HashSet<NormalOrderedSequence<TIndex>>();

            TmpTerms.Push(new LadderSequence<TIndex>(ladderOperator));

            // Anti-commutes creation and annihilation operators to canonical order
            // and creates new terms if spin-orbital indices match.
            while (TmpTerms.Any())
            {
                var tmpTerm = TmpTerms.Pop();
                if (tmpTerm.IsInNormalOrder())
                {
                    NewTerms.Add(new NormalOrderedSequence<TIndex>(tmpTerm));
                }
                else
                {
                    // Anticommute creation and annihilation operators.
                    for (int i = 0; i < tmpTerm.Sequence.Count() - 1; i++)
                    {
                        if ((int)tmpTerm.Sequence.ElementAt(i).Type > (int)tmpTerm.Sequence.ElementAt(i + 1).Type)
                        {
                            // If the two elements have the same spin orbital index, generate a new term.
                            if (tmpTerm.Sequence.ElementAt(i).Index.Equals(tmpTerm.Sequence.ElementAt(i + 1).Index))
                            {
                                var newTerm = new LadderSequence<TIndex>(tmpTerm);
                                newTerm.Sequence.RemoveRange(i, 2);
                                TmpTerms.Push(newTerm);
                            }

                            // Swap the two elements and flip sign of the coefficient.
                            var tmpOp = tmpTerm.Sequence.ElementAt(i + 1);
                            tmpTerm.Sequence[i + 1] = tmpTerm.Sequence.ElementAt(i);
                            tmpTerm.Sequence[i] = tmpOp;
                            tmpTerm.Coefficient *= -1;
                        }
                    }
                    TmpTerms.Push(tmpTerm);
                }
            }
            return NewTerms;
        }

        /// <summary>
        ///  Converts a <see cref="LadderSequence"/> to normal order, then index order. 
        ///  In general, this can generate new terms and modifies the coefficient.
        /// </summary>
        public static HashSet<IndexOrderedSequence<TIndex>> ToIndexOrder<TIndex>(this LadderSequence<TIndex> ladderOperator)
            where TIndex : IEquatable<TIndex>, IComparable<TIndex>
        {
            return new HashSet<IndexOrderedSequence<TIndex>>(
                ladderOperator.ToNormalOrder().Select(o => new IndexOrderedSequence<TIndex>(o))
                );
        }
        #endregion
    }

}



