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
        public List<LadderOperator> Sequence { get; set; } = new List<LadderOperator>();

        /// <summary>
        /// sign (-1,+1) coefficient of ladder operators.
        /// </summary>
        public int Coefficient { get; set; } = 0;

        #region Constructors
        /// <summary>
        /// Constructor for empty ladder operator sequence.
        /// </summary>B
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
            Sequence = ladderOperators.Sequence.ToList();
            Coefficient = ladderOperators.Coefficient;
        }

        /// <summary>
        /// Construct instance from sequence of ladder operators.
        /// </summary>
        /// <param name="setSequence">Sequence of ladder operators.</param>
        public LadderSequence(IEnumerable<LadderOperator> setSequence, int setCoefficient = 1)
        {
            Sequence = setSequence.ToList();
            Coefficient = setCoefficient;
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
            return Sequence.Count() == 0 ? true : Sequence.Select(o => (int)o.Type).IsIntArrayAscending();
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
            return new LadderSequence(left.Sequence.Concat(right.Sequence), left.Coefficient * right.Coefficient);
        }

        /// <summary>
        /// Anti-commutation of ladder operators {x,y}.
        /// </summary>
        /// <param name="x">Left ladder operator.</param>
        /// <param name="y">Right ladder operator.</param>
        /// <returns>Result of {x,y}.</returns>
        public virtual (RaisingLowering, int) AntiCommutator(LadderOperator x, LadderOperator y)
        {
            // {a_x, a_y^\dag} = \delta_{xy}
            // {a_x, a_y} = 1
            // {a_x^\dag, a_y^\dag} = 1
            if (x.Index == y.Index)
            {
                if (x.Type != y.Type)
                {
                    if (x.Type == RaisingLowering.d)
                    {
                        return (RaisingLowering.identity, 1);
                    }
                    else
                    {
                        return (RaisingLowering.identity, -1);
                    }
                }
                else
                {
                    return (RaisingLowering.identity, 0);
                }
            }
            return (RaisingLowering.identity, 0);
        }

        /// <summary>
        /// Counts the number of unique system indices across all <see cref="LadderOperator"/> terms
        /// in a <see cref="LadderSequence"/>
        /// </summary>
        /// <returns>Number of unique system indices.</returns>
        public int GetUniqueIndices()
        {
            return Sequence.Select(o => o.Index).Distinct().Count();
        }

        /// <summary>
        /// Returns a copy of the ladder sequence base class.
        /// </summary>
        /// <returns>LadderSequence class of underlying sequence of ladder operators.</returns>
        public LadderSequence GetLadderSequence()
        {
            return new LadderSequence(Sequence);
        }

        /// <summary>
        /// Returns list of indices of ladder operator sequence.
        /// </summary>
        /// <returns>LadderSequence class of underlying sequence of ladder operators.</returns>
        public List<int> GetLadderSequenceIndices()
        {
            return Sequence.Select(o => o.Index).ToList();
        }

        /// <summary>
        /// Returns list of raising and lowering types of ladder operator sequence.
        /// </summary>
        /// <returns>LadderSequence class of underlying sequence of ladder operators.</returns>
        public List<RaisingLowering> GetLadderSequenceRaisingLowering()
        {
            return Sequence.Select(o => o.Type).ToList();
        }

        /// <summary>
        /// Returns a human-readable description this object.
        /// </summary>
        public override string ToString() 
        {
            return $"{Coefficient} * {string.Join(" ",Sequence)}";
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
            return Sequence.SequenceEqual(x.Sequence);
        }

        public override int GetHashCode()
        {
            int h = 19;
            foreach (var i in Sequence)
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



