// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using System.Linq;
using System.Collections.Generic;

namespace Microsoft.Quantum.Chemistry
{
    /// <summary>
    /// Class representing a sequence of fermionic raising and lowering operators, subject to the additional constraints: 
    /// 1) Normal-ordered, where all raising operators are to the left of all lowering operators.
    /// 2) Index-ordered, where are raising(lowering) operators are in ascending(descending) order.
    /// </summary>
    public class FermionTerm : IndexOrderedLadderSequence
    {
        #region Constructors
        /// <summary>
        /// Constructor for empty instance.
        /// </summary>
        internal FermionTerm() : base() { }

        /// <summary>
        /// Construct a copy of the input instance.
        /// </summary>
        /// <param name="term">Sequence of ladder operators.</param>
        internal FermionTerm(FermionTerm term)
        {
            // All constructions are pass by value.
            sequence = term.sequence.ToList();
            coefficient = term.coefficient;
        }

        /// <summary>
        /// Construct instance from a normal-ordered sequence of ladder operators.
        /// </summary>
        /// <param name="ladderOperators">Normal-ordered sequence of ladder operators.</param>
        public FermionTerm(LadderSequence ladderOperators) : base(ladderOperators) { }
        #endregion

        // If LadderSequence is used as a base class for non-fermionic operators in the future,
        // will need to override the anticommutator.
        //public override (LadderOperator.Type, int) AntiCommutator(LadderOperator x, LadderOperator y)
        //{
        //    return base.AntiCommutator(x, y);
        //}
    }
}



