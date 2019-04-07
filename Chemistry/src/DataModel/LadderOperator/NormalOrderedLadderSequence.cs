// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using System;
using System.Linq;

namespace Microsoft.Quantum.Chemistry.LadderOperators
{
    /// <summary>
    /// Class representing a sequence of raising and lowering operators, subject to the additional constraints: 
    /// 1) Normal-ordered, where all raising operators are to the left of all lowering operators.
    /// </summary>
    public class NormalOrderedLadderSequence : LadderSequence
    {
        #region Constructors
        /// <summary>
        /// Constructor for empty ladder operator sequence.
        /// </summary>
        internal NormalOrderedLadderSequence() : base() { }

        /// <summary>
        /// Construct a copy of the input instance.
        /// </summary>
        /// <param name="ladderOperators">Sequence of ladder operators.</param>
        public NormalOrderedLadderSequence(NormalOrderedLadderSequence ladderOperators)
        {
            // All constructions are pass by value.
            sequence = ladderOperators.sequence.Select(o => o).ToList();
            coefficient = ladderOperators.coefficient;
        }

        /// <summary>
        /// Construct instance from a normal-ordered sequence of ladder operators.
        /// </summary>
        /// <param name="ladderOperators">Normal-ordered sequence of ladder operators.</param>
        public NormalOrderedLadderSequence(LadderSequence ladderOperators) : base(ladderOperators)
        {
            ExceptionIfNotInNormalOrder();
        }
        #endregion

        #region Ordering testers
        /// <summary>
        ///  Checks if raising operators indices are in ascending order, 
        ///  then if lowering operator indices are in descending order.
        /// </summary>
        /// <returns>
        /// Returns <c>true</c> this condition is satisfied.
        /// Returns <c>false</c> otherwise.
        /// </returns>
        public bool IsInIndexOrder()
        {
            return IsInIndexCreationCanonicalOrder() && IsInIndexAnnihilationCanonicalOrder();
        }

        /// <summary>
        ///  Checks whether the creation operator sequence of a <see cref="LadderSequence"/> is in 
        ///  canonical order. This means
        ///  1) <c>SpinOrbital</c> is sorted in ascending order for the creation operators.
        /// </summary>
        /// <returns><c>true</c> if the creation opeartor sequence of a <c>LadderSequence</c> is in 
        /// canonical order. <c>false</c> otherwise</returns>
        public bool IsInIndexCreationCanonicalOrder()
        {
            return sequence.Where(o => o.type == LadderType.u).Select(o => o.index).IsIntArrayAscending();
        }

        /// <summary>
        ///  Checks whether the annihilation operator sequence of a <see cref="LadderSequence"/> is in 
        ///  canonical order. This means
        ///  1) <c>SpinOrbital</c> is sorted in descending order for the annihilation operators.
        /// </summary>
        /// <returns><c>true</c> if the annihilation opeartor sequence of a <c>LadderSequence</c> is in 
        /// canonical order. <c>false</c> otherwise</returns>
        public bool IsInIndexAnnihilationCanonicalOrder()
        {
            return sequence.Where(o => o.type == LadderType.d).Select(o => o.index).Reverse().IsIntArrayAscending();
        }
        #endregion

        #region Reordering methods
        /// <summary>
        ///  Converts a <see cref="NormalOrderedLadderSequence"/> to index order. 
        ///  In general, this can generate new terms and modifies the coefficient.
        /// </summary>
        public void ToIndexOrder()
        {
            var tmp = new NormalOrderedLadderSequence(this);
            if (!tmp.IsInIndexOrder())
            {
                var left = tmp.sequence.Select((op, idx) => new { op, idx }).Where(x => x.op.type == LadderType.u);
                var right = tmp.sequence.Select((op, idx) => new { op, idx }).Where(x => x.op.type == LadderType.d);

                var upArrayIndices = tmp.sequence.Select((op, idx) => new { op, idx }).Where(x => x.op.type == LadderType.u).Select(x => x.idx).ToArray();
                var downArrayIndices = tmp.sequence.Select((op, idx) => new { op, idx }).Where(x => x.op.type == LadderType.d).Select(x => x.idx).ToArray();

                // Bubble sort spin-orbital indices of creation operator.
                while (!tmp.IsInIndexCreationCanonicalOrder())
                {
                    for (int idx = 0; idx < upArrayIndices.Count() - 1; idx++)
                    {
                        if (tmp.sequence.ElementAt(upArrayIndices.ElementAt(idx)).index > tmp.sequence.ElementAt(upArrayIndices.ElementAt(idx + 1)).index)
                        {
                            var tmpLadderOperator = tmp.sequence.ElementAt(upArrayIndices.ElementAt(idx));
                            tmp.sequence[upArrayIndices.ElementAt(idx)] = tmp.sequence[upArrayIndices.ElementAt(idx + 1)];
                            tmp.sequence[upArrayIndices.ElementAt(idx + 1)] = tmpLadderOperator;
                            tmp.coefficient = -1 * tmp.coefficient;
                        }
                    }
                }

                // Bubble sort spin-orbital indices of annihilation operator.
                while (!tmp.IsInIndexAnnihilationCanonicalOrder())
                {
                    for (int idx = 0; idx < downArrayIndices.Length - 1; idx++)
                    {
                        if (tmp.sequence.ElementAt(downArrayIndices.ElementAt(idx)).index < tmp.sequence.ElementAt(downArrayIndices.ElementAt(idx + 1)).index)
                        {
                            var tmpLadderOperator = tmp.sequence.ElementAt(downArrayIndices.ElementAt(idx));
                            tmp.sequence[downArrayIndices.ElementAt(idx)] = tmp.sequence[downArrayIndices.ElementAt(idx + 1)];
                            tmp.sequence[downArrayIndices.ElementAt(idx + 1)] = tmpLadderOperator;
                            tmp.coefficient = -1 * tmp.coefficient;
                        }
                    }
                }
            }
            sequence = tmp.sequence;
            coefficient = tmp.coefficient;
        }
        #endregion

        /// <summary>
        /// This throws an ArgumentException if the operators in NormalOrderedLadderSequence are not normal-ordered.
        /// </summary>
        private void ExceptionIfNotInNormalOrder()
        {
            if (!base.IsInNormalOrder())
            {
                throw new ArgumentException("NormalOrderedLadderSequence must contatin normal-ordered LadderSequence");
            }
        }
    }
    
}



