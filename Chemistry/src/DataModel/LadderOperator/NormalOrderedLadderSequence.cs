// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using System;
using System.Linq;
using System.Collections.Generic;

namespace Microsoft.Quantum.Chemistry
{
    
    public class NormalOrderedLadderSequence : LadderSequence
    {
        #region Constructors
        /// <summary>
        /// Constructor for empty ladder operator sequence.
        /// </summary>
        internal NormalOrderedLadderSequence() : base() { }

        /// <summary>
        /// Construct <see cref="NormalOrderedLadderSequence"/> from another <see cref="NormalOrderedLadderSequence"/>.
        /// </summary>
        /// <param name="setSequence">Sequence of ladder operators.</param>
        public NormalOrderedLadderSequence(NormalOrderedLadderSequence ladderOperators)
        {
            // All constructions are pass by value.
            sequence = ladderOperators.sequence.Select(o => o).ToList();
            coefficient = ladderOperators.coefficient;
        }

        /// <summary>
        /// Construct <see cref="NormalOrderedLadderSequence"/> from sequence of ladder operators.
        /// </summary>
        /// <param name="setSequence">Sequence of ladder operators.</param>
        public NormalOrderedLadderSequence(LadderSequence ladderOperators) : base(ladderOperators)
        {
            ExceptionIfNotInNormalOrder();
        }

        /// <summary>
        /// Construct <see cref="NormalOrderedLadderSequence"/> from sequence of ladder operators.
        /// </summary>
        /// <param name="setSequence">Sequence of ladder operators.</param>
        public NormalOrderedLadderSequence(IEnumerable<LadderOperator> setSequence, int setCoefficient = 1) : base(setSequence, setCoefficient)
        {
            ExceptionIfNotInNormalOrder();
        }


        /// <summary>
        /// Construct <see cref="NormalOrderedLadderSequence"/> from sequence of ladder operators.
        /// </summary>
        /// <param name="set">Sequence of ladder operators.</param>
        public NormalOrderedLadderSequence(IEnumerable<(LadderOperator.Type, int)> set) : base(set)
        {
            ExceptionIfNotInNormalOrder();
        }

        /// <summary>
        /// FermionTerm constructor that assumes normal-ordered fermionic 
        /// creation and annihilation operators, and that the number of
        /// creation an annihilation operators are equal.
        /// </summary>
        public NormalOrderedLadderSequence(IEnumerable<int> indices) : base(indices)
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
            return sequence.Where(o => o.type == LadderOperator.Type.u).Select(o => o.index).IsIntArrayAscending();
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
            return sequence.Where(o => o.type == LadderOperator.Type.d).Select(o => o.index).Reverse().IsIntArrayAscending();
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



