// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using System;
using System.Linq;
using System.Collections.Generic;

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
    public class NormalOrderedSequence<TIndex> : LadderSequence<TIndex>
        where TIndex : IEquatable<TIndex>
    {
        #region Constructors
        /// <summary>
        /// Constructor for empty ladder operator sequence.
        /// </summary>
        internal NormalOrderedSequence() : base() { }

        /// <summary>
        /// Construct instance from a normal-ordered sequence of ladder operators.
        /// </summary>
        /// <param name="ladderOperators">Normal-ordered sequence of ladder operators.</param>
        public NormalOrderedSequence(LadderSequence<TIndex> ladderOperators) : base(ladderOperators) => ThrowExceptionIfNotInNormalOrder();

        /// <summary>
        /// Construct instance from a normal-ordered sequence of ladder operators.
        /// </summary>
        /// <param name="ladderOperators">Normal-ordered sequence of ladder operators.</param>
        public NormalOrderedSequence(IEnumerable<LadderOperator<TIndex>> ladderOperators, int setSign = 1) : base(ladderOperators, setSign) => ThrowExceptionIfNotInNormalOrder();
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
        public bool IsInIndexCreationCanonicalOrder() => Sequence.Where(o => o.GetRaisingLowering() == RaisingLowering.u).Select(o => o.GetIndex()).IsInAscendingOrder();

        /// <summary>
        ///  Checks whether the annihilation operator sequence of a <see cref="LadderSequence"/> is in 
        ///  canonical order. This means
        ///  1) <c>SpinOrbital</c> is sorted in descending order for the annihilation operators.
        /// </summary>
        /// <returns><c>true</c> if the annihilation opeartor sequence of a <c>LadderSequence</c> is in 
        /// canonical order. <c>false</c> otherwise</returns>
        public bool IsInIndexAnnihilationCanonicalOrder() => Sequence.Where(o => o.GetRaisingLowering() == RaisingLowering.d).Select(o => o.GetIndex()).Reverse().IsInAscendingOrder();
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



