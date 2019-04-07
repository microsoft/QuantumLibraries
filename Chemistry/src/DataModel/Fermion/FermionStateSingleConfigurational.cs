// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using System;
using System.Linq;
using System.Collections.Generic;

using Microsoft.Quantum.Chemistry.LadderOperators;

namespace Microsoft.Quantum.Chemistry.Fermion
{
    // This class is for Fermion terms that are not grouped into Hermitian bunches.
    // Maybe need a stype for quantum state?

    /// <summary>
    /// Class representing a sequence of fermionic raising and lowering operators, subject to the additional constraints: 
    /// 1) Normal-ordered, where all raising operators are to the left of all lowering operators.
    /// 2) Index-ordered, where are raising(lowering) operators are in ascending(descending) order.
    /// 3) Contains only creation operators.
    /// </summary>
    public class FermionStateSingleConfigurational : FermionTerm
    {
        #region Constructors
        /// <summary>
        /// Constructor for empty instance.
        /// </summary>
        internal FermionStateSingleConfigurational() : base() { }

        /// <summary>
        /// Construct a copy of the input instance.
        /// </summary>
        /// <param name="term">Sequence of ladder operators.</param>
        internal FermionStateSingleConfigurational(FermionStateSingleConfigurational term)
        {
            // All constructions are pass by value.
            sequence = term.sequence.Select(o => o).ToList();
            coefficient = term.coefficient;
        }

        /// <summary>
        /// Construct instance from a normal-ordered sequence of ladder operators.
        /// </summary>
        /// <param name="ladderOperators">Normal-ordered sequence of ladder operators.</param>
        public FermionStateSingleConfigurational(LadderSequence ladderOperators) : base(ladderOperators) { }
        #endregion

        /// <summary>
        /// This throws an ArgumentException if the operators in NormalOrderedLadderSequence are not normal-ordered.
        /// </summary>
        private void ExceptionIfNotOnlyRaising()
        {
            if (sequence.Where(o => o.type == LadderType.d).Count() > 0)
            {
                throw new ArgumentException("FermionStateSingleConfigurational must contatin only raising operators.");
            }
        }

    }

}



