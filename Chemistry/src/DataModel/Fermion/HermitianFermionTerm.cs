// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using System;
using System.Linq;
using System.Collections.Generic;

namespace Microsoft.Quantum.Chemistry
{
    public class HermitianFermionTerm : IndexOrderedLadderOperators, HamiltonianTerm<TermType.Fermion>
    {
        internal HermitianFermionTerm() : base() { }

        // Disallow inputs that are not normal ordered.
        internal HermitianFermionTerm(HermitianFermionTerm term)
        {
            sequence = term.sequence;
            coefficient = term.coefficient;
        }

        public HermitianFermionTerm(IEnumerable<LadderOperator> setSequence, int setCoefficient = 1) : base(setSequence, setCoefficient) { ToCanonicalOrder(); }
        public HermitianFermionTerm(IEnumerable<(LadderOperator.Type, int)> set) : base(set) { ToCanonicalOrder(); }
        public HermitianFermionTerm(LadderOperatorSequence set) : base(set) { ToCanonicalOrder(); }
        public HermitianFermionTerm(IEnumerable<int> indices) : base(indices) { ToCanonicalOrder(); }

        /// <summary>
        ///  Checks if raising operators indices are in ascending order, 
        ///  then if lowering operator indices are in descending order.
        /// </summary>
        /// <returns>
        /// Returns <c>true</c> this condition is satisfied.
        /// Returns <c>false</c> otherwise.
        /// </returns>
        public bool IsInCanonicalOrder()
        {
            if (base.IsInIndexOrder())
            {
                var creationSequence = sequence.Where(o => o.type == LadderOperator.Type.u).Select(o => o.index);
                var annihilationSequence = sequence.Where(o => o.type == LadderOperator.Type.d).Select(o => o.index);
                if (creationSequence.Count() == annihilationSequence.Count())
                {
                    if (Extensions.CompareIntArray(creationSequence, annihilationSequence.Reverse()) > 0)
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

        private void ToCanonicalOrder()
        {
            // Take Hermitian Conjugate    
            if (!IsInCanonicalOrder())
            {
                sequence = sequence.Select(o => (o.type == LadderOperator.Type.d ? LadderOperator.Type.u : LadderOperator.Type.d, o.index)).Select(o => new LadderOperator(o)).Reverse().ToList();
            }
        }

        public TermType.Fermion GetTermType()
        {
            var length = sequence.Count();
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
        
    }
    
}



