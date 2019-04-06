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
    public class FermionTermSingle : IndexOrderedLadderSequence
    {
        #region Constructors
        /// <summary>
        /// Constructor for empty ladder operator sequence.
        /// </summary>
        internal FermionTermSingle() : base() { }

        /// <summary>
        /// Construct a copy of the input instance.
        /// </summary>
        /// <param name="term">Fermion term.</param>
        internal FermionTermSingle(FermionTermSingle term)
        {
            // All constructions are pass by value.
            sequence = term.sequence.Select(o => o).ToList();
            coefficient = term.coefficient;
        }

        public FermionTermSingle(LadderSequence set) : base(set) { }
        #endregion
    }
}



