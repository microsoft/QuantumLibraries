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
        public FermionTerm(List<LadderOperator> setSequence, Symmetry setSymmetry) : base(setSequence) { symmetry = setSymmetry; }
        public FermionTerm(List<(LadderOperator.Type, int)> set, Symmetry setSymmetry) : base(set) { }
        public FermionTerm(LadderOperators set, Symmetry setSymmetry) : this(set.sequence, setSymmetry) { }


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
        private bool IsInCanonicalOrder()
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


        public (int, FermionTerm) ToCanonicalOrderFromNormalOrder(bool AllowHermitianConjugate = false)
        {
            var (sign, newTerm) = base.ToCanonicalOrderFromNormalOrder();

            // Take Hermitian conjugate if still not in canonical order. 
            if (!newTerm.IsInIndexOrder())
            {
                newTerm.sequence = newTerm.sequence.Select(o => (o.type == LadderOperator.Type.d ? LadderOperator.Type.u : LadderOperator.Type.d, o.index)).Select(o => new LadderOperator(o)).Reverse().ToList();
            }
            return (sign, new FermionTerm(newTerm));
        }

        /// <summary>
        ///  Converts a <c>FermionTerm</c> to canonical order. This generates
        ///  new terms and modifies the coefficient as needed.
        /// </summary>
        public List<(int, FermionTerm)> ToCanonicalOrder()
        {
            // Step 1: anti-commute creation to the left.
            // Step 2: sort to canonical order

            var TmpTerms = new Stack<(int, FermionTerm)>();
            var NewTerms = new List<(int, FermionTerm)>();

            TmpTerms.Push((1,this));

            // Anti-commutes creation and annihilation operators to canonical order
            // and creates new terms if spin-orbital indices match.
            while (TmpTerms.Any())
            {
                var tmpTerm = TmpTerms.Pop();
                if (tmpTerm.Item2.IsInNormalOrder())
                {
                    NewTerms.Add(tmpTerm);
                }
                else
                {
                    // Anticommute creation and annihilation operators.
                    for (int i = 0; i < tmpTerm.Item2.sequence.Count() - 1; i++)
                    {
                        if ((int) tmpTerm.Item2.sequence.ElementAt(i).type > (int) tmpTerm.Item2.sequence.ElementAt(i + 1).type)
                        {
                            // Swap the two elements and flip sign of the coefficient.
                            tmpTerm.Item2.sequence[i + 1] = tmpTerm.Item2.sequence.ElementAt(i);
                            tmpTerm.Item2.sequence[i] = tmpTerm.Item2.sequence.ElementAt(i + 1);

                            var antiCommutedTerm = (-1 * tmpTerm.Item1, new FermionTerm(tmpTerm.Item2));

                            TmpTerms.Push(antiCommutedTerm);

                            // If the two elements have the same spin orbital index, generate a new term.
                            if (tmpTerm.Item2.sequence.ElementAt(i).type == tmpTerm.Item2.sequence.ElementAt(i + 1).type)
                            {
                                tmpTerm.Item2.sequence.RemoveRange(i, 2);

                                var newTerm = (tmpTerm.Item1, new FermionTerm(tmpTerm.Item2));

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
                var (coeff, tmp) = NewTerms[idx];
                var (coeff2, tmp2) = tmp.ToCanonicalOrderFromNormalOrder();
                NewTerms[idx] = (coeff * coeff2, tmp2);
            }
            return NewTerms;
        }

    }
}



