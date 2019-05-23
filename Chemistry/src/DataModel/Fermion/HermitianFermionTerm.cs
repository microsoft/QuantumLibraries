// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using System;
using System.Linq;
using System.Collections.Generic;

using Microsoft.Quantum.Chemistry.LadderOperators;
using Microsoft.Quantum.Chemistry;
using Microsoft.Quantum.Chemistry;

namespace Microsoft.Quantum.Chemistry.Fermion
{
    using static Microsoft.Quantum.Chemistry.Extensions;
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
    /// <item>
    /// <description>
    /// Hermitian, and is assumed to be implicitly summed with its Hermitian conjugate if not explicitly Hermitian.
    /// </description>
    /// </item>
    /// </summary>
    public class AntiHermitianFermionTerm : HermitianFermionTerm
    {
        #region Constructors
        /// <summary>
        /// Constructor for empty instance.
        /// </summary>
        public AntiHermitianFermionTerm() : base() { }
        #endregion

        /// <summary>
        /// Additional sort using the Hermitian conjugate.
        /// </summary>
        private void NormalizeToCanonicalOrder()
        {
            // Take Hermitian Conjugate    
            if (!IsInCanonicalOrder())
            {
                Sequence = Sequence.Select(o => (o.Type == RaisingLowering.d ? RaisingLowering.u : RaisingLowering.d, o.Index)).Select(o => new FermionOperator(o)).Reverse().ToList();
                Coefficient *= -1;
            }
        }
    }


}



