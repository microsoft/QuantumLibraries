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

    public struct LadderOperator
    {
        public enum Type
        {
            u = 0, d = 1;
        }

        public Type type;
        public int index;

        public int Commutation(LadderOperator x, LadderOperator y)
        {
            if(x.index == y.index)
            {
                if(x.type != y.type)
                {
                    if(x.type == Type.d){
                        return 1;
                    }
                    else
                    {
                        return -1;
                    }
                }
            }
            return 0;
        }

        public LadderOperator(Type setType, int setIndex)
        {
            type = setType;
            this.index = setIndex;
        }

        public LadderOperator((Type, int) set) : this(set.Item1, set.Item2) { }

        public override bool Equals(object obj)
        {
            return (obj is LadderOperator x) ? Equals(x) : false;
        }

        public bool Equals(LadderOperator x)
        {
            return (type == x.type) && (index == x.index);
        }

        public override int GetHashCode()
        {
            return index * Enum.GetNames(typeof(Type)).Length + (int) type;
        }

        public static bool operator ==(LadderOperator x, LadderOperator y)
        {
            return x.Equals(y);
        }

        public static bool operator !=(LadderOperator x, LadderOperator y)
        {
            return !(x.Equals(y));
        }
    }

    public struct LadderOperators
    {
        public IEnumerable<LadderOperator> sequence;

        public LadderOperators(IEnumerable<LadderOperator> setSequence)
        {
            sequence = setSequence;
        }

        public LadderOperators(IEnumerable<(LadderOperator.Type, int)> set) : this(set.Select(o => new LadderOperator(o))) { }

        public bool IsInNormalOrder(LadderOperators x)
        {
            return x.sequence.Count() == 0 ? true : x.sequence.Select(o => (int) o.type).IsIntArrayAscending();
        }

        public IEnumerable<(int, LadderOperators)> NormalOrdering(LadderOperators x)
        {


            // Step 1: anti-commute creation to the left.
            // Step 2: sort to canonical order

            var TmpTerms = new Stack<(int,LadderOperators)>();
            var NewTerms = new List<(int, LadderOperators)>();

            TmpTerms.Push(x);

            // Anti-commutes creation and annihilation operators to canonical order
            // and creates new terms if spin-orbital indices match.
            while (TmpTerms.Any())
            {
                var tmpTerm = TmpTerms.Pop();
                if (tmpTerm.IsInNormalOrder())
                {
                    NewTerms.Add(tmpTerm);
                }
                else
                {
                    // Anticommute creation and annihilation operators.
                    for (int i = 0; i < tmpTerm.CreationAnnihilationIndices.Count() - 1; i++)
                    {
                        if (tmpTerm.CreationAnnihilationIndices.ElementAt(i) < tmpTerm.CreationAnnihilationIndices.ElementAt(i + 1))
                        {
                            var antiCommutedCreationAnnihilationIndices = tmpTerm.CreationAnnihilationIndices.ToList();
                            var antiCommutedSpinOrbitalIndices = tmpTerm.SpinOrbitalIndices.ToList();
                            // Swap the two elements and flip sign of the coefficient.
                            antiCommutedCreationAnnihilationIndices[i + 1] = tmpTerm.CreationAnnihilationIndices.ElementAt(i);
                            antiCommutedCreationAnnihilationIndices[i] = tmpTerm.CreationAnnihilationIndices.ElementAt(i + 1);
                            antiCommutedSpinOrbitalIndices[i + 1] = tmpTerm.SpinOrbitalIndices.ElementAt(i);
                            antiCommutedSpinOrbitalIndices[i] = tmpTerm.SpinOrbitalIndices.ElementAt(i + 1);

                            var antiCommutedTerm = new FermionTerm(antiCommutedCreationAnnihilationIndices.ToArray(), antiCommutedSpinOrbitalIndices.ToArray(), -1.0 * tmpTerm.coeff);

                            TmpTerms.Push(antiCommutedTerm);

                            // If the two elements have the same spin orbital index, generate a new term.
                            if (antiCommutedSpinOrbitalIndices.ElementAt(i).Equals(antiCommutedSpinOrbitalIndices.ElementAt(i + 1)))
                            {
                                var newCreationAnnihilationIndices = antiCommutedCreationAnnihilationIndices.ToList();
                                var newSpinOrbitalIndices = antiCommutedSpinOrbitalIndices.ToList();
                                newCreationAnnihilationIndices.RemoveRange(i, 2);
                                newSpinOrbitalIndices.RemoveRange(i, 2);

                                var newTerm = new FermionTerm(newCreationAnnihilationIndices.ToArray(), newSpinOrbitalIndices.ToArray(), tmpTerm.coeff);

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
                var tmp = NewTerms[idx];
                tmp.ToSpinOrbitalCanonicalOrder();
                NewTerms[idx] = tmp;
            }
            return NewTerms;
        }

        public override bool Equals(object obj)
        {
            return (obj is LadderOperators x) ? Equals(x) : false;
        }

        public bool Equals(LadderOperators x)
        {
            return sequence.SequenceEqual(x.sequence);
        }

        public override int GetHashCode()
        {
            int h = 19;
            foreach (var i in sequence)
            {
                h = h * 53 + i.GetHashCode();
            }
            return h;
        }

        public static bool operator ==(LadderOperators x, LadderOperators y)
        {
            return x.Equals(y);
        }

        public static bool operator !=(LadderOperators x, LadderOperators y)
        {
            return !(x.Equals(y));
        }
    }
}



