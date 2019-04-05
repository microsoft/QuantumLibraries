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
            u = 0, d = 1
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

    public partial class FermionTerm : IEquatable<FermionTerm>
    {
        public List<LadderOperator> sequence { get; set; }

        public FermionTerm(List<LadderOperator> setSequence)
        {
            sequence = setSequence;
        }

        public FermionTerm(List<(LadderOperator.Type, int)> set) : this(set.Select(o => new LadderOperator(o)).ToList()) { }

        public bool IsInNormalOrder()
        {
            return sequence.Count() == 0 ? true : sequence.Select(o => (int) o.type).IsIntArrayAscending();
        }

        /// <summary>
        ///  Checks whether the creation operator sequence of a <see cref="FermionTerm"/> is in 
        ///  canonical order. This means
        ///  1) <c>SpinOrbital</c> is sorted in ascending order for the creation operators.
        /// </summary>
        /// <returns><c>true</c> if the creation opeartor sequence of a <c>FermionTerm</c> is in 
        /// canonical order. <c>false</c> otherwise</returns>
        private bool IsInIndexCreationCanonicalOrder()
        {
            return sequence.Where(o => o.type == LadderOperator.Type.u).Select(o => o.index).IsIntArrayAscending();
        }

        /// <summary>
        ///  Checks whether the annihilation operator sequence of a <see cref="FermionTerm"/> is in 
        ///  canonical order. This means
        ///  1) <c>SpinOrbital</c> is sorted in descending order for the annihilation operators.
        /// </summary>
        /// <returns><c>true</c> if the annihilation opeartor sequence of a <c>FermionTerm</c> is in 
        /// canonical order. <c>false</c> otherwise</returns>
        private bool IsInIndexAnnihilationCanonicalOrder()
        {
            return sequence.Where(o => o.type == LadderOperator.Type.d).Select(o => o.index).Reverse().IsIntArrayAscending(); 
        }

        public bool IsInCanonicalOrder()
        {
            var creation = sequence.Where(o => o.type == LadderOperator.Type.u).Select(o => (long) o.index);
            var annihilation = sequence.Where(o => o.type == LadderOperator.Type.d).Select(o => (long) o.index);
            if (creation.Count() == annihilation.Count())
            {
                if (Extensions.CompareIntArray(creation, annihilation.Reverse()) > 0)
                {
                    return false;
                }
            }
            return IsInNormalOrder()
                && IsInIndexCreationCanonicalOrder()
                && IsInIndexAnnihilationCanonicalOrder();

        }

        public (int, FermionTerm) ToCanonicalOrder(bool AllowHermitianConjugate = true)
        {
            FermionTerm tmp = new FermionTerm(sequence);
            int coeff = 1;
            //var left = sequence.Where(o => o.type == LadderOperator.Type.u).OrderBy(o => o.index);
            //var right = sequence.Where(o => o.type == LadderOperator.Type.d).OrderBy(o => o.index).Reverse();
            //var normalOrdered = new LadderOperators(left.Concat(right));
            //return (0, normalOrdered);

            // Check that FermionTerm is normal-ordered.
            if (!tmp.IsInNormalOrder())
            {
                throw new System.ArgumentException(
                    $"ToCanonicalOrder() assumes input is normal-ordered. This is currently not satisfied."
                    );
            }
            // Check that FermionTerm spin-orbital indices are in canonical order.
            if (!tmp.IsInCanonicalOrder())
            {

                var upArrayIndices = tmp.sequence.Select((op, idx) => new { op, idx }).Where(x => x.op.type == LadderOperator.Type.u).Select(x => x.idx).ToArray();
                var downArrayIndices = tmp.sequence.Select((op, idx) => new { op, idx }).Where(x => x.op.type == LadderOperator.Type.d).Select(x => x.idx).ToArray();

                // Bubble sort spin-orbital indices of creation operator.
                while (!tmp.IsInIndexCreationCanonicalOrder())
                {
                    for (int idx = 0; idx < upArrayIndices.Count() - 1; idx++)
                    {
                        if (tmp.sequence.ElementAt(upArrayIndices.ElementAt(idx)).index > tmp.sequence.ElementAt(upArrayIndices.ElementAt(idx + 1)).index)
                        {
                            var tmpLadderOperator = tmp.sequence.ElementAt(upArrayIndices.ElementAt(idx));
                            tmp.sequence[upArrayIndices.ElementAt(idx)] = tmp.sequence[upArrayIndices.ElementAt(idx + 1)];
                            tmp.sequence[upArrayIndices.ElementAt(idx + 1)] = tmpLadderOperator;
                            coeff = -1 * coeff;
                        }
                    }
                }

                // Bubble sort spin-orbital indices of annihilation operator.
                while (!tmp.IsInIndexAnnihilationCanonicalOrder())
                {
                    for (int idx = 0; idx < downArrayIndices.Length - 1; idx++)
                    {
                        if (tmp.sequence.ElementAt(downArrayIndices.ElementAt(idx)).index < tmp.sequence.ElementAt(downArrayIndices.ElementAt(idx + 1)).index)
                        {
                            var tmpLadderOperator = tmp.sequence.ElementAt(downArrayIndices.ElementAt(idx));
                            tmp.sequence[downArrayIndices.ElementAt(idx)] = tmp.sequence[downArrayIndices.ElementAt(idx + 1)];
                            tmp.sequence[downArrayIndices.ElementAt(idx + 1)] = tmpLadderOperator;
                            coeff = -1 * coeff;
                        }
                    }
                }
            }

            // Take Hermitian conjugate if still not in canonical order. 
            if (!tmp.IsInCanonicalOrder())
            {
                tmp.sequence = tmp.sequence.Select(o => (o.type == LadderOperator.Type.d ? LadderOperator.Type.u : LadderOperator.Type.d, o.index)).Select(o => new LadderOperator(o)).Reverse().ToList();
            }
            return (coeff, tmp);
        }

        public override bool Equals(object obj)
        {
            return (obj is FermionTerm x) ? Equals(x) : false;
        }

        public bool Equals(FermionTerm x)
        {
            // If parameter is null, return false.
            if (ReferenceEquals(x, null))
            {
                return false;
            }

            // Optimization for a common success case.
            if (ReferenceEquals(this, x))
            {
                return true;
            }

            // If run-time types are not exactly the same, return false.
            if (GetType() != x.GetType())
            {
                return false;
            }
            // Return true if the fields match.
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

        public static bool operator == (FermionTerm x, FermionTerm y)
        {
            // Check for null on left side.
            if (Object.ReferenceEquals(x, null))
            {
                if (Object.ReferenceEquals(y, null))
                {
                    // null == null = true.
                    return true;
                }

                // Only the left side is null.
                return false;
            }
            // Equals handles case of null on right side.
            return x.Equals(y);
        }

        public static bool operator !=(FermionTerm x, FermionTerm y)
        {
            return !(x == y);
        }
    }   

    public partial class FermionTerm
    {
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


    }


}



