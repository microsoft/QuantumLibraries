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
        public FermionTerm(LadderOperators set) : this(set.sequence) { }
        public FermionTerm(List<LadderOperator> setSequence) : base(setSequence) { }
        public FermionTerm(List<(LadderOperator.Type, int)> set) : base() { }

        /// <summary>
        /// FermionTerm constructor that assumes normal-ordered fermionic 
        /// creation and annihilation operators, and that the number of
        /// creation an annihilation operators are equal.
        /// </summary>
        public FermionTerm(IEnumerable<int> indices, bool sort = true)
        {
            var length = indices.Count();
            if (length % 2 == 1)
            {
                throw new System.ArgumentException(
                    $"SpinOrbital array of length {length} must be of even length."
                    );
            }
            var tmp = new FermionTerm(indices.Select((o, idx) => new LadderOperator((idx < length / 2 ? LadderOperator.Type.u : LadderOperator.Type.d, o))).ToList());

            sequence = tmp.sequence;
        }

        public Int64 GetUniqueIndices()
        {
            return sequence.Select(o => o.index).Distinct().Count();
        }

        public (int, FermionTerm) ToCanonicalOrderFromNormalOrder(bool AllowHermitianConjugate = true)
        {
            var (coeff, tmp) = base.ToCanonicalOrderFromNormalOrder();

            // Take Hermitian conjugate if still not in canonical order. 
            if (!tmp.IsInCanonicalOrder())
            {
                tmp.sequence = tmp.sequence.Select(o => (o.type == LadderOperator.Type.d ? LadderOperator.Type.u : LadderOperator.Type.d, o.index)).Select(o => new LadderOperator(o)).Reverse().ToList();
            }
            return (coeff, new FermionTerm(tmp));
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



