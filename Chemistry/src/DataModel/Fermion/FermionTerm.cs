// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using System.Linq;
using System;
using System.Collections.Generic;
using System.Linq;

using Microsoft.Quantum.Chemistry.Generic;

namespace Microsoft.Quantum.Chemistry.Fermion
{
    using FermionOperator = LadderOperator<int>;

    /// <summary>
    /// Class representing a sequence of fermionic raising and lowering operators, subject to the additional constraints: 
    /// <list type="number">
    /// <item>
    /// <description>
    /// Normal-ordered, where all raising operators are to the left of all lowering operators.
    /// </description>
    /// </item>
    /// <item>
    /// <description>
    /// Index-ordered, where are raising(lowering) operators are in ascending(descending) order.
    /// </description>
    /// </item>
    /// </summary>
    public class FermionTerm : IndexOrderedSequence<int>
    {
        #region Constructors
        /// <summary>
        /// Constructor for empty instance.
        /// </summary>
        internal FermionTerm() : base() { }

        /// <summary>
        /// Construct fermion term instance from a normal-ordered sequence of ladder operators.
        /// </summary>
        /// <param name="ladderOperators">Normal-ordered sequence of ladder operators.</param>
        public FermionTerm(LadderSequence<int> ladderOperators) : base(ladderOperators) { }

        /// <summary>
        /// Construct instance from a normal-ordered sequence of ladder operators.
        /// </summary>
        /// <param name="ladderOperators">Normal-ordered sequence of ladder operators.</param>
        public FermionTerm(IEnumerable<FermionOperator> ladderOperators, int setSign = 1 ) : base(ladderOperators, setSign) { }
        #endregion

        // This exists as a convenience function for creating fermion terms in samples.
        /// <summary>
        /// Implicit operator for creating a Ladder operator.
        /// </summary>
        /// <param name="setOperator">Tuple where the first parameter
        /// is the raising or lowering index, and the second parameter
        /// is the position index of the ladder operator.</param>
        public static implicit operator FermionTerm((RaisingLowering, int)[] setSequence)
        {
            return new FermionTerm(setSequence.Select(o => new FermionOperator(o)));
        }

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
        public static implicit operator FermionTerm(int[] indices)
        {
            var length = indices.Count();
            if (length % 2 == 1)
            {
                throw new System.ArgumentException(
                    $"Number of terms provided is `{length}` and must be of even length."
                    );
            }
            Func<int, int, (RaisingLowering, int)> GetLadderOperator = (index, position)
                => (position < length / 2 ? RaisingLowering.u : RaisingLowering.d, index);
            return indices.Select((o, idx) => GetLadderOperator(o, idx)).ToArray();
        }
    }
}



