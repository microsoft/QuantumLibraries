// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using System.Linq;

using Microsoft.Quantum.Chemistry.LadderOperators;

namespace Microsoft.Quantum.Chemistry.Fermion
{

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
            Sequence = term.Sequence.ToList();
            Coefficient = term.Coefficient;
        }

        /// <summary>
        /// Construct instance from a normal-ordered sequence of ladder operators.
        /// </summary>
        /// <param name="ladderOperators">Normal-ordered sequence of ladder operators.</param>
        public FermionTerm(LadderSequence ladderOperators) : base(ladderOperators) { }
        #endregion

        // If LadderSequence is used as a base class for non-fermionic operators in the future,
        // will need to override the anticommutator.
        //public override (LadderType, int) AntiCommutator(LadderOperator x, LadderOperator y)
        //{
        //    return base.AntiCommutator(x, y);
        //}
    }
}



