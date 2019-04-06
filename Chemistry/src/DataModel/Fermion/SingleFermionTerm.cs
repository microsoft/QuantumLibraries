// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using System.Collections.Generic;

namespace Microsoft.Quantum.Chemistry
{
    // This class is for Fermion terms that are not grouped into Hermitian bunches.
    
    public class SingleFermionTerm : IndexOrderedLadderOperators
    {
        internal SingleFermionTerm() : base() { }

        // Disallow inputs that are not normal ordered.
        internal SingleFermionTerm(SingleFermionTerm term)
        {
            sequence = term.sequence;
            coefficient = term.coefficient;
        }

        public SingleFermionTerm(IEnumerable<LadderOperator> setSequence, int setCoefficient = 1) : base(setSequence, setCoefficient) { }
        public SingleFermionTerm(IEnumerable<(LadderOperator.Type, int)> set) : base(set) { }
        public SingleFermionTerm(LadderOperatorSequence set) : base(set) { }
        public SingleFermionTerm(IEnumerable<int> indices) : base(indices) { }
        
    }
}



