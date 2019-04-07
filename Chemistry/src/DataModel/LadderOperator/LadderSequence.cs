// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using System;
using System.Linq;
using System.Collections.Generic;

namespace Microsoft.Quantum.Chemistry.LadderOperators
{
    
    /// <summary>
    /// Class representing a sequence of raising and lowering operators.
    /// </summary>
    public class LadderSequence : IEquatable<LadderSequence>
    {

        /// <summary>
        /// Sequence of ladder operators.
        /// </summary>
        public List<LadderOperator> sequence { get; set; } = new List<LadderOperator>();

        /// <summary>
        /// sign (-1,+1) coefficient of ladder operators.
        /// </summary>
        public int coefficient { get; set; } = 0;

        #region Constructors
        /// <summary>
        /// Constructor for empty ladder operator sequence.
        /// </summary>
        public LadderSequence() {
        }

        /// <summary>
        /// Construct a copy of the input instance.
        /// </summary>
        /// <param name="ladderOperators">Sequence of ladder operators.</param>
        public LadderSequence(LadderSequence ladderOperators)
        {
            // All constructions are pass by value. 
            // ToList() copies values if the underlying object is a value type.
            sequence = ladderOperators.sequence.ToList();
            coefficient = ladderOperators.coefficient;
        }

        /// <summary>
        /// Construct instance from sequence of ladder operators.
        /// </summary>
        /// <param name="setSequence">Sequence of ladder operators.</param>
        public LadderSequence(IEnumerable<LadderOperator> setSequence, int setCoefficient = 1)
        {
            sequence = setSequence.ToList();
            coefficient = setCoefficient;
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
        /// <param name="left">Left <see cref="LadderSequence"/> <c>x</c>.</param>
        /// <param name="right">Right <see cref="LadderSequence"/> <c>y</c>.</param>
        /// <returns>
        /// Returns new <see cref="LadderSequence"/> <c>xy</c> where coefficients and 
        /// LadderOperatorSequences are multipled together.
        /// </returns>
        public LadderSequence Multiply(LadderSequence left, LadderSequence right)
        {
            return new LadderSequence(left.sequence.Concat(right.sequence), left.coefficient * right.coefficient);
        }

        /// <summary>
        /// Anti-commutation of ladder operators {x,y}.
        /// </summary>
        /// <param name="x">Left ladder operator.</param>
        /// <param name="y">Right ladder operator.</param>
        /// <returns>Result of {x,y}.</returns>
        public virtual (LadderType, int) AntiCommutator(LadderOperator x, LadderOperator y)
        {
            // {a_x, a_y^\dag} = \delta_{xy}
            // {a_x, a_y} = 1
            // {a_x^\dag, a_y^\dag} = 1
            if (x.index == y.index)
            {
                if (x.type != y.type)
                {
                    if (x.type == LadderType.d)
                    {
                        return (LadderType.identity, 1);
                    }
                    else
                    {
                        return (LadderType.identity, -1);
                    }
                }
                else
                {
                    return (LadderType.identity, 0);
                }
            }
            return (LadderType.identity, 0);
        }

        /// <summary>
        /// Counts the number of unique system indices across all <see cref="LadderOperator"/> terms
        /// in a <see cref="LadderSequence"/>
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
            return (obj is LadderSequence x) ? Equals(x) : false;
        }

        public bool Equals(LadderSequence x)
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

        public static bool operator == (LadderSequence x, LadderSequence y)
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

        public static bool operator !=(LadderSequence x, LadderSequence y)
        {
            return !(x == y);
        }
        #endregion


    }

}



