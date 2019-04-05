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
            u = 0, d = 1, identity
        }

        
        public Type type;
        public int index;
        

        public static (Type, int) AntiCommutation(LadderOperator x, LadderOperator y)
        {
            // {a_x, a_y^\dag} = \delta_{xy}
            // {a_x, a_y} = 1
            // {a_x^\dag, a_y^\dag} = 1
            if (x.index == y.index)
            {
                if(x.type != y.type)
                {
                    if(x.type == Type.d){
                        return (Type.identity, 1);
                    }
                    else
                    {
                        return (Type.identity, -1);
                    }
                }
                else
                {
                    return (Type.identity, 0);
                }
            }
            return (Type.identity, 0);
        }

        public LadderOperator(Type setType, int setIndex)
        {
            type = setType;
            index = setIndex;
        }

        public LadderOperator((Type, int) set) : this(set.Item1, set.Item2) { }

        #region Equality Testing
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
#endregion
    }


    public class LadderOperators : IEquatable<LadderOperators>
    {
        public enum Ordering
        {
            NoOrder = 0, NormalOrder = 1, CanonicalOrder = 2
        }
        private Ordering Order;
        public Ordering order
        {
            get
            {
                return Order;
            }
            set
            {
                throw new ArgumentException("order may not be changed directly");
            }
        }

        private List<LadderOperator> Sequence;

        public List<LadderOperator> sequence
        {
            get
            {
                return Sequence;
            }
            set
            {
                Sequence = value;
                order = GetOrdering();
            }
        }

        /// <summary>
        /// Ordering of terms. 
        /// </summary>

        

        public LadderOperators()
        {
        }

        public LadderOperators(List<LadderOperator> setSequence)
        {
            sequence = setSequence;
        }

        public LadderOperators(List<(LadderOperator.Type, int)> set) : this(set.Select(o => new LadderOperator(o)).ToList()) { }

        /// <summary>
        /// FermionTerm constructor that assumes normal-ordered fermionic 
        /// creation and annihilation operators, and that the number of
        /// creation an annihilation operators are equal.
        /// </summary>
        public LadderOperators(IEnumerable<int> indices)
        {
            var length = indices.Count();
            if (length % 2 == 1)
            {
                throw new System.ArgumentException(
                    $"Number of terms provided is `{length}` and must be of even length."
                    );
            }
            var tmp = new LadderOperators(indices.Select((o, idx) => new LadderOperator((idx < length / 2 ? LadderOperator.Type.u : LadderOperator.Type.d, o))).ToList());

            sequence = tmp.sequence;
        }

        /// <summary>
        /// Checks whether all raising operators are to the left of all lowering operators.
        /// </summary>
        /// <returns>
        /// Returns <c>true</c> this condition is satisfied.
        /// Returns <c>false</c> otherwise.
        /// </returns>
        private bool IsInNormalOrder()
        {
            bool isInNormalOrder = sequence.Count() == 0 ? true : sequence.Select(o => (int)o.type).IsIntArrayAscending();
            return isInNormalOrder;
        }

        /// <summary>
        ///  Checks if raising operators indices are in ascending order, 
        ///  then if lowering operator indices are in descending order.
        /// </summary>
        /// <returns>
        /// Returns <c>true</c> this condition is satisfied.
        /// Returns <c>false</c> otherwise.
        /// </returns>
        private bool IsInIndexOrder()
        {
            return IsInIndexCreationCanonicalOrder() && IsInIndexAnnihilationCanonicalOrder();
        }

        public void SetOrdering()
        {
            Order = GetOrdering();
        }

        /// <summary>
        ///  Checks whether the ladder operators are normal ordered, then checks if raising operators indices
        ///  are in ascending order, then if lowering operator indices are in descending order.
        /// </summary>
        /// <returns>
        /// Returns the ordering type detected.
        /// </returns>
        public Ordering GetOrdering()
        { 
            if (!IsInNormalOrder())
            {
                return Ordering.NoOrder;
            }
            else if (!IsInIndexOrder())
            {
                return Ordering.NormalOrder;
            }
            else
            {
                return Ordering.CanonicalOrder;
            }
        }

        /// <summary>
        ///  Checks whether the creation operator sequence of a <see cref="LadderOperators"/> is in 
        ///  canonical order. This means
        ///  1) <c>SpinOrbital</c> is sorted in ascending order for the creation operators.
        /// </summary>
        /// <returns><c>true</c> if the creation opeartor sequence of a <c>LadderOperators</c> is in 
        /// canonical order. <c>false</c> otherwise</returns>
        private bool IsInIndexCreationCanonicalOrder()
        {
            return sequence.Where(o => o.type == LadderOperator.Type.u).Select(o => o.index).IsIntArrayAscending();
        }

        /// <summary>
        ///  Checks whether the annihilation operator sequence of a <see cref="LadderOperators"/> is in 
        ///  canonical order. This means
        ///  1) <c>SpinOrbital</c> is sorted in descending order for the annihilation operators.
        /// </summary>
        /// <returns><c>true</c> if the annihilation opeartor sequence of a <c>LadderOperators</c> is in 
        /// canonical order. <c>false</c> otherwise</returns>
        private bool IsInIndexAnnihilationCanonicalOrder()
        {
            return sequence.Where(o => o.type == LadderOperator.Type.d).Select(o => o.index).Reverse().IsIntArrayAscending(); 
        }

        public Int64 GetUniqueIndices()
        {
            return sequence.Select(o => o.index).Distinct().Count();
        }



        public (LadderOperators, int) ToCanonicalOrderFromNormalOrder()
        {
            var tmp = new LadderOperators(sequence);
            int coeff = 1;
            //var left = sequence.Where(o => o.type == LadderOperator.Type.u).OrderBy(o => o.index);
            //var right = sequence.Where(o => o.type == LadderOperator.Type.d).OrderBy(o => o.index).Reverse();
            //var normalOrdered = new LadderOperators(left.Concat(right));
            //return (0, normalOrdered);

            // Check that LadderOperators is normal-ordered.
            if (!tmp.IsInNormalOrder())
            {
                throw new System.ArgumentException(
                    $"ToCanonicalOrder() assumes input is normal-ordered. This is currently not satisfied."
                    );
            }
            // Check that LadderOperators spin-orbital indices are in canonical order.
            if (!tmp.IsInIndexOrder())
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
            tmp = new LadderOperators(tmp.Sequence);
            return (tmp, coeff);
        }


        #region Equality Testing

        public override bool Equals(object obj)
        {
            return (obj is LadderOperators x) ? Equals(x) : false;
        }

        public bool Equals(LadderOperators x)
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

        public static bool operator == (LadderOperators x, LadderOperators y)
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

        public static bool operator !=(LadderOperators x, LadderOperators y)
        {
            return !(x == y);
        }
        #endregion
    }
}



