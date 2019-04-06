// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using System.Linq;
using System.Collections.Generic;

namespace Microsoft.Quantum.Chemistry
{
    // This class is for Fermion terms that are not grouped into Hermitian bunches.
    
    public class SingleFermionTerm : IndexOrderedLadderSequence
    {
        internal SingleFermionTerm() : base() { }

        // Disallow inputs that are not normal ordered.
        internal SingleFermionTerm(SingleFermionTerm term)
        {
            // All constructions are pass by value.
            sequence = term.sequence.Select(o => o).ToList();
            coefficient = term.coefficient;
        }

        public SingleFermionTerm(IEnumerable<LadderOperator> setSequence, int setCoefficient = 1) : base(setSequence, setCoefficient) { }
        public SingleFermionTerm(IEnumerable<(LadderOperator.Type, int)> set) : base(set) { }
        public SingleFermionTerm(LadderSequence set) : base(set) { }
        public SingleFermionTerm(IEnumerable<int> indices) : base(indices) { }
        
    }
}



