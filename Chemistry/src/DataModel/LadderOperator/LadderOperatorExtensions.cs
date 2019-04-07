// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using System;
using System.Linq;
using System.Collections.Generic;

namespace Microsoft.Quantum.Chemistry.LadderOperators
{
    using static LadderType;

    public static partial class Extensions
    {
        #region Convenience constructors
        /// <summary>
        /// Construct <see cref="LadderSequence"/> from sequence of ladder operators.
        /// </summary>
        /// <param name="setSequence">Sequence of ladder operators.</param>
        /// <param name="setSign">Set the sign coefficient of the sequence.</param>
        /// <returns>
        /// Sequence of ladder operators.
        /// </returns>
        public static LadderSequence ToLadderSequence(this IEnumerable<(LadderType, int)> setSequence, int setSign = 1) {
            return new LadderSequence(setSequence.Select(o => new LadderOperator(o)), setSign);
            var tmp = new[] { (u, 1), (u, 2), (d, 3), (d, 4) }.ToLadderSequence();
        }

        /// <summary>
        /// Construct <see cref="LadderSequence"/> from an even-length sequence of integers.
        /// </summary>
        /// <param name="setSequence">Even-length sequence of integers.</param>
        /// <param name="setSign">Set the sign coefficient of the sequence.</param>
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
        public static LadderSequence ToLadderSequence(this IEnumerable<int> indices)
        {
            var length = indices.Count();
            if (length % 2 == 1)
            {
                throw new System.ArgumentException(
                    $"Number of terms provided is `{length}` and must be of even length."
                    );
            }
            Func<int, int, (LadderType, int)> GetLadderOperator = (index, position)
                => (position < length / 2 ? u : d, index);
            return indices.Select((o, idx) => GetLadderOperator(o, idx)).ToLadderSequence();
        }
        #endregion

        #region Reordering methods

        /// <summary>
        ///  Converts a <see cref="LadderSequence"/> to normal order. 
        ///  In general, this can generate new terms and modifies the coefficient.
        /// </summary>
        public static HashSet<NormalOrderedLadderSequence> CreateNormalOrder(this LadderSequence ladderOperator)
        {
            // Recursively anti-commute creation to the left.
            var TmpTerms = new Stack<LadderSequence>();
            var NewTerms = new HashSet<NormalOrderedLadderSequence>();

            TmpTerms.Push(new LadderSequence(ladderOperator));

            // Anti-commutes creation and annihilation operators to canonical order
            // and creates new terms if spin-orbital indices match.
            while (TmpTerms.Any())
            {
                var tmpTerm = TmpTerms.Pop();
                if (tmpTerm.IsInNormalOrder())
                {
                    NewTerms.Add(new NormalOrderedLadderSequence(tmpTerm));
                }
                else
                {
                    // Anticommute creation and annihilation operators.
                    for (int i = 0; i < tmpTerm.sequence.Count() - 1; i++)
                    {
                        if ((int)tmpTerm.sequence.ElementAt(i).type > (int)tmpTerm.sequence.ElementAt(i + 1).type)
                        {
                            // If the two elements have the same spin orbital index, generate a new term.
                            if (tmpTerm.sequence.ElementAt(i).index == tmpTerm.sequence.ElementAt(i + 1).index)
                            {
                                var newTerm = new LadderSequence(tmpTerm);
                                newTerm.sequence.RemoveRange(i, 2);
                                TmpTerms.Push(newTerm);
                            }

                            // Swap the two elements and flip sign of the coefficient.
                            var tmpOp = tmpTerm.sequence.ElementAt(i + 1);
                            tmpTerm.sequence[i + 1] = tmpTerm.sequence.ElementAt(i);
                            tmpTerm.sequence[i] = tmpOp;
                            tmpTerm.coefficient *= -1;
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
        public static HashSet<IndexOrderedLadderSequence> CreateIndexOrder(this LadderSequence ladderOperator)
        {
            return new HashSet<IndexOrderedLadderSequence>(
                ladderOperator.CreateNormalOrder().Select(o => new IndexOrderedLadderSequence(o))
                );
        }
        #endregion
    }

}



