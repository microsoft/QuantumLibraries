// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.


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
    /// <item>
    /// <description>
    /// Index-ordered, where are raising(lowering) operators are in ascending(descending) order.
    /// </description>
    /// </item>
    /// </summary>
    public class IndexOrderedLadderSequence : NormalOrderedLadderSequence
    {
        #region Constructors
        /// <summary>
        /// Constructor for empty ladder operator sequence.
        /// </summary>
        internal IndexOrderedLadderSequence() : base() { }

        /// <summary>
        /// Construct a copy of the input instance.
        /// </summary>
        /// <param name="ladderOperators">Sequence of ladder operators.</param>
        public IndexOrderedLadderSequence(IndexOrderedLadderSequence ladderOperators)
        {
            // All constructions are pass by value.
            Sequence = ladderOperators.Sequence.ToList();
            Coefficient = ladderOperators.Coefficient;
        }

        /// <summary>
        /// Construct instance from a normal-ordered sequence of ladder operators.
        /// </summary>
        /// <param name="ladderOperators">Normal-ordered sequence of ladder operators.</param>
        public IndexOrderedLadderSequence(LadderSequence ladderOperators) : base(ladderOperators)
        {
            NormalizeToIndexOrder();
        }
        #endregion
        
    }

}



