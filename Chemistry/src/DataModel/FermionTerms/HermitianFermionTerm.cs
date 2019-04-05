// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using Microsoft.Quantum.Simulation.Core;

using System;
using System.Linq;
using System.Collections.Generic;
using Microsoft.Quantum.Chemistry;
using System.Numerics;

namespace Microsoft.Quantum.Chemistry
{
    public class HermitianFermionTerm : NormalOrderedLadderOperators
    {
        internal HermitianFermionTerm() : base() { }

        // Disallow inputs that are not normal ordered.
        public HermitianFermionTerm(HermitianFermionTerm term)
        {
            sequence = term.sequence;
            coefficient = term.coefficient;
        }

        public HermitianFermionTerm(List<LadderOperator> setSequence, int setCoefficient = 1) : base(setSequence, setCoefficient) { }
        public HermitianFermionTerm(List<(LadderOperator.Type, int)> set) : base(set) { }
        public HermitianFermionTerm(LadderOperators set) : base(set) { }
        public HermitianFermionTerm(IEnumerable<int> indices) : base(indices) { }

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
                return true;
            }
            return false;
        }


        public HermitianFermionTerm CreateCanonicalOrderOrder()
        {
            var fermionTerm = new HermitianFermionTerm(base.CreateCreateIndexOrder());

            // Take Hermitian Conjugate
            if (!fermionTerm.IsInCanonicalOrder())
            {
                fermionTerm.sequence = fermionTerm.sequence.Select(o => (o.type == LadderOperator.Type.d ? LadderOperator.Type.u : LadderOperator.Type.d, o.index)).Select(o => new LadderOperator(o)).Reverse().ToList();
            }
            return fermionTerm;
        }
    }
    
}



