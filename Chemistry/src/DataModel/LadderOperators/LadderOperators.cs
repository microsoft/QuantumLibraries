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
    
    public class LadderOperators : IEquatable<LadderOperators>
    {


        public List<LadderOperator> sequence { get; set; }

        public int coefficient { get; set; }

        #region Constructors
        internal LadderOperators() { }

        /// <summary>
        /// Construct LadderOperators from another LadderOperators.
        /// </summary>
        /// <param name="setSequence">Sequence of ladder operators.</param>
        public LadderOperators(LadderOperators ladderOperators)
        {
            sequence = ladderOperators.sequence;
            coefficient = ladderOperators.coefficient;
        }

        /// <summary>
        /// Construct LadderOperators from sequence of ladder operators.
        /// </summary>
        /// <param name="setSequence">Sequence of ladder operators.</param>
        public LadderOperators(IEnumerable<LadderOperator> setSequence, int setCoefficient = 1)
        {
            sequence = setSequence.ToList();
            coefficient = setCoefficient;
        }

        /// <summary>
        /// Construct LadderOperators from sequence of ladder operators.
        /// </summary>
        /// <param name="setSequence">Sequence of ladder operators.</param>
        public LadderOperators(IEnumerable<(LadderOperator.Type, int)> set, int setCoefficient = 1) : this(set.Select(o => new LadderOperator(o)).ToList(), setCoefficient) { }
        

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
            Func<int, int, LadderOperator> GetLadderOperator = (index, position) 
                => new LadderOperator(position < length / 2 ? LadderOperator.Type.u : LadderOperator.Type.d, index);
            var tmp = new LadderOperators(indices.Select((o, idx) => GetLadderOperator(o,idx)).ToList());

            sequence = tmp.sequence;
        }
        #endregion
        
        /// <summary>
        /// Checks whether all raising operators are to the left of all lowering operators.
        /// </summary>
        /// <returns>
        /// Returns <c>true</c> this condition is satisfied.
        /// Returns <c>false</c> otherwise.
        /// </returns>
        public bool IsInNormalOrder()
        {
            return sequence.Count() == 0 ? true : sequence.Select(o => (int)o.type).IsIntArrayAscending();
        }


        public Int64 GetUniqueIndices()
        {
            return sequence.Select(o => o.index).Distinct().Count();
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



