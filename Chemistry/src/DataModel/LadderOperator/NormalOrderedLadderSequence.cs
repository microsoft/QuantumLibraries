// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using System;
using System.Linq;

namespace Microsoft.Quantum.Chemistry.LadderOperators
{
    /// <summary>
    /// Class representing a sequence of raising and lowering operators, subject to the additional constraints: 
    /// <list type="number">
    /// <item>
    /// <description>
    /// Normal-ordered, where all raising operators are to the left of all lowering operators.
    /// </description>
    /// </item>
    /// </item>
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
            if (ladderOperators != null)
            {
                // All constructions are pass by value.
                Sequence = ladderOperators.Sequence.ToList();
                Coefficient = ladderOperators.Coefficient;
            }
        }

        /// <summary>
        /// Construct instance from a normal-ordered sequence of ladder operators.
        /// </summary>
        /// <param name="ladderOperators">Normal-ordered sequence of ladder operators.</param>
        public NormalOrderedLadderSequence(LadderSequence ladderOperators) : base(ladderOperators) => ThrowExceptionIfNotInNormalOrder();
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
        public bool IsInIndexOrder() => IsInIndexCreationCanonicalOrder() && IsInIndexAnnihilationCanonicalOrder();

        /// <summary>
        ///  Checks whether the creation operator sequence of a <see cref="LadderSequence"/> is in 
        ///  canonical order. This means
        ///  1) <c>SpinOrbital</c> is sorted in ascending order for the creation operators.
        /// </summary>
        /// <returns><c>true</c> if the creation opeartor sequence of a <c>LadderSequence</c> is in 
        /// canonical order. <c>false</c> otherwise</returns>
        public bool IsInIndexCreationCanonicalOrder() => Sequence.Where(o => o.Type == RaisingLowering.u).Select(o => o.Index).IsInAscendingOrder();

        /// <summary>
        ///  Checks whether the annihilation operator sequence of a <see cref="LadderSequence"/> is in 
        ///  canonical order. This means
        ///  1) <c>SpinOrbital</c> is sorted in descending order for the annihilation operators.
        /// </summary>
        /// <returns><c>true</c> if the annihilation opeartor sequence of a <c>LadderSequence</c> is in 
        /// canonical order. <c>false</c> otherwise</returns>
        public bool IsInIndexAnnihilationCanonicalOrder() => Sequence.Where(o => o.Type == RaisingLowering.d).Select(o => o.Index).Reverse().IsInAscendingOrder();
        #endregion

        #region Reordering methods
        /// <summary>
        ///  Converts a <see cref="NormalOrderedLadderSequence"/> to index order. 
        ///  In general, this can generate new terms and modifies the coefficient.
        /// </summary>
        public void NormalizeToIndexOrder()
        {
            var tmp = new NormalOrderedLadderSequence(this);
            if (!tmp.IsInIndexOrder())
            {
                var left = tmp.Sequence.Select((op, idx) => new { op, idx }).Where(x => x.op.Type == RaisingLowering.u);
                var right = tmp.Sequence.Select((op, idx) => new { op, idx }).Where(x => x.op.Type == RaisingLowering.d);

                var upArrayIndices = tmp.Sequence.Select((op, idx) => new { op, idx }).Where(x => x.op.Type == RaisingLowering.u).Select(x => x.idx).ToArray();
                var downArrayIndices = tmp.Sequence.Select((op, idx) => new { op, idx }).Where(x => x.op.Type == RaisingLowering.d).Select(x => x.idx).ToArray();

                // Bubble sort spin-orbital indices of creation operator.
                while (!tmp.IsInIndexCreationCanonicalOrder())
                {
                    for (int idx = 0; idx < upArrayIndices.Count() - 1; idx++)
                    {
                        if (tmp.Sequence.ElementAt(upArrayIndices.ElementAt(idx)).Index > tmp.Sequence.ElementAt(upArrayIndices.ElementAt(idx + 1)).Index)
                        {
                            var tmpLadderOperator = tmp.Sequence.ElementAt(upArrayIndices.ElementAt(idx));
                            tmp.Sequence[upArrayIndices.ElementAt(idx)] = tmp.Sequence[upArrayIndices.ElementAt(idx + 1)];
                            tmp.Sequence[upArrayIndices.ElementAt(idx + 1)] = tmpLadderOperator;
                            tmp.Coefficient = -1 * tmp.Coefficient;
                        }
                    }
                }

                // Bubble sort spin-orbital indices of annihilation operator.
                while (!tmp.IsInIndexAnnihilationCanonicalOrder())
                {
                    for (int idx = 0; idx < downArrayIndices.Length - 1; idx++)
                    {
                        if (tmp.Sequence.ElementAt(downArrayIndices.ElementAt(idx)).Index < tmp.Sequence.ElementAt(downArrayIndices.ElementAt(idx + 1)).Index)
                        {
                            var tmpLadderOperator = tmp.Sequence.ElementAt(downArrayIndices.ElementAt(idx));
                            tmp.Sequence[downArrayIndices.ElementAt(idx)] = tmp.Sequence[downArrayIndices.ElementAt(idx + 1)];
                            tmp.Sequence[downArrayIndices.ElementAt(idx + 1)] = tmpLadderOperator;
                            tmp.Coefficient = -1 * tmp.Coefficient;
                        }
                    }
                }
            }
            Sequence = tmp.Sequence;
            Coefficient = tmp.Coefficient;
        }
        #endregion

        /// <summary>
        /// This throws an ArgumentException if the operators in NormalOrderedLadderSequence are not normal-ordered.
        /// </summary>
        private void ThrowExceptionIfNotInNormalOrder()
        {
            if (!base.IsInNormalOrder())
            {
                throw new ArgumentException("NormalOrderedLadderSequence must contatin normal-ordered LadderSequence");
            }
        }
    }

}



