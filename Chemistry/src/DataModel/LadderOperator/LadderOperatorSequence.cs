// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using System;
using System.Linq;
using System.Collections.Generic;

namespace Microsoft.Quantum.Chemistry
{
    
    public class LadderOperatorSequence : IEquatable<LadderOperatorSequence>
    {

        /// <summary>
        /// Sequence of ladder operators.
        /// </summary>
        public List<LadderOperator> sequence { get; set; }

        /// <summary>
        /// sign (-1,+1) coefficient of ladder operators.
        /// </summary>
        public int coefficient { get; set; }

        #region Constructors
        /// <summary>
        /// Constructor for empty ladder operator sequence.
        /// </summary>
        internal LadderOperatorSequence() { }

        /// <summary>
        /// Construct LadderOperators from another LadderOperators.
        /// </summary>
        /// <param name="setSequence">Sequence of ladder operators.</param>
        public LadderOperatorSequence(LadderOperatorSequence ladderOperators)
        {
            // All constructions are pass by value.
            sequence = ladderOperators.sequence.Select(o => o).ToList();
            coefficient = ladderOperators.coefficient;
        }

        /// <summary>
        /// Construct LadderOperators from sequence of ladder operators.
        /// </summary>
        /// <param name="setSequence">Sequence of ladder operators.</param>
        public LadderOperatorSequence(IEnumerable<LadderOperator> setSequence, int setCoefficient = 1)
        {
            sequence = setSequence.Select(o => o).ToList();
            coefficient = setCoefficient;
        }

        /// <summary>
        /// Construct LadderOperators from sequence of ladder operators.
        /// </summary>
        /// <param name="setSequence">Sequence of ladder operators.</param>
        public LadderOperatorSequence(IEnumerable<(LadderOperator.Type, int)> set, int setCoefficient = 1) : this(set.Select(o => new LadderOperator(o)).ToList(), setCoefficient) { }
        

        /// <summary>
        /// FermionTerm constructor that assumes normal-ordered fermionic 
        /// creation and annihilation operators, and that the number of
        /// creation an annihilation operators are equal.
        /// </summary>
        public LadderOperatorSequence(IEnumerable<int> indices)
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
            var tmp = new LadderOperatorSequence(indices.Select((o, idx) => GetLadderOperator(o,idx)).ToList());

            sequence = tmp.sequence;
            coefficient = tmp.coefficient;
        }
        #endregion

        #region Ordering testers
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
        #endregion

        /// <summary>
        /// Concatenates two Fermion terms.
        /// </summary>
        /// <param name="left">Left <see cref="LadderOperatorSequence"/> <c>x</c>.</param>
        /// <param name="right">Right <see cref="LadderOperatorSequence"/> <c>y</c>.</param>
        /// <returns>
        /// Returns new <see cref="LadderOperatorSequence"/> <c>xy</c> where coefficients and 
        /// LadderOperatorSequences are multipled together.
        /// </returns>
        public LadderOperatorSequence Multiply(LadderOperatorSequence left, LadderOperatorSequence right)
        {
            return new LadderOperatorSequence(left.sequence.Concat(right.sequence), left.coefficient * right.coefficient);
        }

        /// <summary>
        /// Counts the number of unique system indices across all <see cref="LadderOperator"/> terms
        /// in a <see cref="LadderOperatorSequence"/>
        /// </summary>
        /// <returns>Number of unique system indices.</returns>
        public int GetUniqueIndices()
        {
            return sequence.Select(o => o.index).Distinct().Count();
        }

        /// <summary>
        /// Returns a human-readable description this object.
        /// </summary>
        public override string ToString() 
        {
            return $"{coefficient} * {string.Join(" ",sequence)}";
        }     

        #region Equality Testing

        public override bool Equals(object obj)
        {
            return (obj is LadderOperatorSequence x) ? Equals(x) : false;
        }

        public bool Equals(LadderOperatorSequence x)
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

        public static bool operator == (LadderOperatorSequence x, LadderOperatorSequence y)
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

        public static bool operator !=(LadderOperatorSequence x, LadderOperatorSequence y)
        {
            return !(x == y);
        }
        #endregion


    }

}



