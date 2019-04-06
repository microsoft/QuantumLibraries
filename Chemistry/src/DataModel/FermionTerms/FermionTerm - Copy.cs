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
    public partial class FermionTerm : LadderOperators
    {
        // FermionTerms that are in normal order should always be sorted into index order.
        // Disallow inputs that are not normal ordered.

        public FermionTerm(List<LadderOperator> setSequence, Symmetry setSymmetry, int setCoefficient = 1) : base(setSequence, setCoefficient) {
            symmetry = setSymmetry;

        }
        public FermionTerm(List<(LadderOperator.Type, int)> set, Symmetry setSymmetry) : base(set) { }
        public FermionTerm(LadderOperators set, Symmetry setSymmetry) : base(set.sequence, set.coefficient) { symmetry = setSymmetry; }
        
        public Symmetry symmetry;
        
        public enum Symmetry
        {
            Single, Hermitian, AntiHermitian 
        }
        
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
            if (symmetry == Symmetry.Single) {
                return base.GetOrdering() == Ordering.CanonicalOrder;
            }
            else if(symmetry == Symmetry.Hermitian)
            {
                Ordering ladderOperatorOrdering = base.GetOrdering();
                if(ladderOperatorOrdering == Ordering.CanonicalOrder)
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
            }
            return false;
        }
        
        /// <summary>
        /// FermionTerm constructor that assumes normal-ordered fermionic 
        /// creation and annihilation operators, and that the number of
        /// creation an annihilation operators are equal.
        /// </summary>
        public FermionTerm(IEnumerable<int> indices, Symmetry setSymmetry) : base(indices)
        {
            symmetry = setSymmetry;
            if(symmetry == Symmetry.Hermitian)
            {

            }
        }


        public FermionTerm CreateCanonicalOrderFromNormalOrder()
        {
            var fermionTerm = new FermionTerm(base.CreateCanonicalOrderFromNormalOrder(), this.symmetry);
            
            if (fermionTerm.symmetry == Symmetry.Hermitian)
            {
                // Take Hermitian Conjugate
                if (!fermionTerm.IsInCanonicalOrder())
                {
                    fermionTerm.sequence = fermionTerm.sequence.Select(o => (o.type == LadderOperator.Type.d ? LadderOperator.Type.u : LadderOperator.Type.d, o.index)).Select(o => new LadderOperator(o)).Reverse().ToList();
                }
            }
            return fermionTerm;
        }

        /// <summary>
        ///  Converts a <c>FermionTerm</c> to canonical order. This generates
        ///  new terms and modifies the coefficient as needed.
        /// </summary>
        public List<FermionTerm> CreateCanonicalOrder()
        {
            // Step 1: anti-commute creation to the left.
            // Step 2: sort to canonical order

            var TmpTerms = new Stack<FermionTerm>();
            var NewTerms = new List<FermionTerm>();

            TmpTerms.Push(this);

            // Anti-commutes creation and annihilation operators to canonical order
            // and creates new terms if spin-orbital indices match.
            while (TmpTerms.Any())
            {
                var tmpTerm = TmpTerms.Pop();
                if (tmpTerm.IsInCanonicalOrder())
                {
                    NewTerms.Add(tmpTerm);
                }
                else
                {
                    // Anticommute creation and annihilation operators.
                    for (int i = 0; i < tmpTerm.sequence.Count() - 1; i++)
                    {
                        if ((int) tmpTerm.sequence.ElementAt(i).type > (int) tmpTerm.sequence.ElementAt(i + 1).type)
                        {
                            // Swap the two elements and flip sign of the coefficient.
                            tmpTerm.sequence[i + 1] = tmpTerm.sequence.ElementAt(i);
                            tmpTerm.sequence[i] = tmpTerm.sequence.ElementAt(i + 1);
                            tmpTerm.coefficient *= -1;
                            var antiCommutedTerm = new FermionTerm(tmpTerm, this.symmetry);

                            TmpTerms.Push(antiCommutedTerm);

                            // If the two elements have the same spin orbital index, generate a new term.
                            if (tmpTerm.sequence.ElementAt(i).type == tmpTerm.sequence.ElementAt(i + 1).type)
                            {
                                tmpTerm.sequence.RemoveRange(i, 2);

                                var newTerm = new FermionTerm(tmpTerm, this.symmetry);

                                TmpTerms.Push(newTerm);
                            }
                            break;
                        }
                    }
                }
            }

            // Anti-commutes spin-orbital indices to canonical order
            // and changes the sign of the coefficient as necessay.
            for (int idx = 0; idx < NewTerms.Count(); idx++)
            {
                NewTerms[idx] = NewTerms[idx].CreateCanonicalOrderFromNormalOrder();
            }
            return NewTerms;
        }

    }
}



