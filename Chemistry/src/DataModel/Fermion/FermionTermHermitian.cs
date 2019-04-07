// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using System;
using System.Linq;
using System.Collections.Generic;

using Microsoft.Quantum.Chemistry.LadderOperators;
using Microsoft.Quantum.Chemistry.Generic;
using Microsoft.Quantum.Chemistry;

namespace Microsoft.Quantum.Chemistry.Fermion
{
    using static Microsoft.Quantum.Chemistry.Extensions;

    /// <summary>
    /// Class representing a sequence of fermionic raising and lowering operators, subject to the additional constraints: 
    /// 1) Normal-ordered, where all raising operators are to the left of all lowering operators.
    /// 2) Index-ordered, where are raising(lowering) operators are in ascending(descending) order.
    /// 3) Hermitian, and is assumed to be implicitly summed with its Hermitian conjugate if not explicitly Hermitian.
    /// </summary>
    public class FermionTermHermitian : FermionTerm, ITermIndex<TermType.Fermion>, IEquatable<FermionTermHermitian>
    {
        #region Constructors
        /// <summary>
        /// Constructor for empty instance.
        /// </summary>
        internal FermionTermHermitian() : base() { }

        /// <summary>
        /// Construct a copy of the input instance.
        /// </summary>
        /// <param name="term">Sequence of ladder operators.</param>
        internal FermionTermHermitian(FermionTermHermitian term)
        {
            // All constructions are pass by value.
            Sequence = term.Sequence.Select(o => o).ToList();
            Coefficient = term.Coefficient;
        }

        /// <summary>
        /// Construct instance from a normal-ordered sequence of ladder operators.
        /// </summary>
        /// <param name="ladderOperators">Normal-ordered sequence of ladder operators.</param>
        public FermionTermHermitian(LadderSequence ladderOperators) : base(ladderOperators) { ToCanonicalOrder(); }
        #endregion

        /// <summary>
        ///  Checks if raising operators indices are in ascending order, 
        ///  then if lowering operator indices are in descending order.
        /// </summary>
        /// <returns>
        /// Returns <c>true</c> this condition is satisfied.
        /// Returns <c>false</c> otherwise.
        /// </returns>
        /// <remarks>
        /// This should always return <c>true</c> when invoked outside this class.
        /// </remarks>
        public bool IsInCanonicalOrder()
        {
            if (IsInIndexOrder())
            {
                var creationSequence = Sequence.Where(o => o.Type == RaisingLowering.u).Select(o => o.Index);
                var annihilationSequence = Sequence.Where(o => o.Type == RaisingLowering.d).Select(o => o.Index);
                if (creationSequence.Count() == annihilationSequence.Count())
                {
                    if (CompareIntArray(creationSequence, annihilationSequence.Reverse()) > 0)
                    {
                        return false;
                    }
                }
                else if(creationSequence.Count() < annihilationSequence.Count())
                {
                    return false;
                }
                return true;
            }
            return false;
        }

        /// <summary>
        /// Additional sort using the Hermitian conjugate.
        /// </summary>
        private void ToCanonicalOrder()
        {
            // Take Hermitian Conjugate    
            if (!IsInCanonicalOrder())
            {
                Sequence = Sequence.Select(o => (o.Type == RaisingLowering.d ? RaisingLowering.u : RaisingLowering.d, o.Index)).Select(o => new LadderOperator(o)).Reverse().ToList();
            }
        }

        /// <summary>
        /// Return the category of this fermion term.
        /// </summary>
        /// <returns>Category of fermion term.</returns>
        public TermType.Fermion GetTermType()
        {
            var length = Sequence.Count();
            var uniqueIndices = this.GetUniqueIndices();

            switch (length)
            {
                case 0:
                    return TermType.Fermion.Identity;
                case 2:
                    switch (uniqueIndices)
                    {
                        case 1:
                            return TermType.Fermion.PP;
                        case 2:
                            return TermType.Fermion.PQ;
                        default:
                            throw new ArgumentException("Attempted to classify unknown fermion term.");
                    }
                case 4:
                    switch (uniqueIndices)
                    {
                        case 2:
                            return TermType.Fermion.PQQP;
                        case 3:
                            return TermType.Fermion.PQQR;
                        case 4:
                            return TermType.Fermion.PQRS;
                        default:
                            throw new ArgumentException("Attempted to classify unknown fermion term.");
                    }
                default:
                    throw new ArgumentException("Attempted to classify unknown fermion term.");
            }
        }

        
        #region Equality Testing

        public override bool Equals(object obj)
        {
            return (obj is FermionTermHermitian x) ? Equals(x) : false;
        }

        public bool Equals(FermionTermHermitian x)
        {
            return base.Equals(x);
        }

        public override int GetHashCode()
        {
            return base.GetHashCode();
        }
        #endregion

    }


}



