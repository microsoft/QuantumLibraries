// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using System;
using System.Linq;
using System.Collections.Generic;

namespace Microsoft.Quantum.Chemistry
{
    
    public class IndexOrderedLadderSequence : NormalOrderedLadderSequence
    {
        #region Constructors
        /// <summary>
        /// Constructor for empty ladder operator sequence.
        /// </summary>
        internal IndexOrderedLadderSequence() : base() { }

        /// <summary>
        /// Construct LadderSequence from another LadderSequence.
        /// </summary>
        /// <param name="setSequence">Sequence of ladder operators.</param>
        public IndexOrderedLadderSequence(IndexOrderedLadderSequence ladderOperators)
        {
            // All constructions are pass by value.
            sequence = ladderOperators.sequence.Select(o => o).ToList();
            coefficient = ladderOperators.coefficient;
        }

        /// <summary>
        /// Construct LadderSequence from sequence of ladder operators.
        /// </summary>
        /// <param name="setSequence">Sequence of ladder operators.</param>
        public IndexOrderedLadderSequence(LadderSequence ladderOperators) : base(ladderOperators)
        {
            ToIndexOrder();
        }
        
        /// <summary>
        /// Construct LadderSequence from sequence of ladder operators.
        /// </summary>
        /// <param name="setSequence">Sequence of ladder operators.</param>
        public IndexOrderedLadderSequence(IEnumerable<LadderOperator> setSequence, int setCoefficient = 1) : base(setSequence, setCoefficient)
        {
            ToIndexOrder();
        }


        /// <summary>
        /// Construct LadderSequence from sequence of ladder operators.
        /// </summary>
        /// <param name="setSequence">Sequence of ladder operators.</param>
        public IndexOrderedLadderSequence(IEnumerable<(LadderOperator.Type, int)> set) : base(set)
        {
            ToIndexOrder();
        }

        /// <summary>
        /// FermionTerm constructor that assumes normal-ordered fermionic 
        /// creation and annihilation operators, and that the number of
        /// creation an annihilation operators are equal.
        /// </summary>
        public IndexOrderedLadderSequence(IEnumerable<int> indices) : base(indices)
        {
            ToIndexOrder();
        }
        #endregion




        #region Reordering methods
        private void ToIndexOrder()
        {
            var ladderTerm = CreateIndexOrder();
            sequence = ladderTerm.sequence;
            coefficient = ladderTerm.coefficient;        
        }

        private IndexOrderedLadderSequence CreateIndexOrder()
        {
            return this.CreateIndexOrder();
        }
        
        #endregion
        
    }

}



