// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using System;
using System.Collections.Generic;
using System.Linq;
using Newtonsoft.Json;
using Microsoft.Quantum.Chemistry;

namespace Microsoft.Quantum.Chemistry.LadderOperators
{


    /// <summary>
    /// Class representing a sequence of raising and lowering operators.
    /// </summary>
    [JsonConverter(typeof(Json.LadderSequenceJsonConverter))]
    public class LadderSequence<TIndex> : IEquatable<LadderSequence<TIndex>>, Json.ILadderSequence
        where TIndex : IEquatable<TIndex>
    {

        /// <summary>
        /// Sequence of ladder operators.
        /// </summary>
        public List<LadderOperator<TIndex>> Sequence { get; set; } = new List<LadderOperator<TIndex>>();

        /// <summary>
        /// sign (-1,+1) coefficient of ladder operators.
        /// </summary>
        public int Coefficient { get; set; } = 1;

        #region Json serialization
        /// <summary>
        /// Returns ladder operator sequence.
        /// </summary>
        /// <returns>Ladder operator sequence.</returns>
        public object _JsonGetSequence() => Sequence;
         /// <summary>
         /// Sets ladder operator sequence.
         /// </summary>
         /// <param name="set">Set ladder operator sequence to this input.</param>
        public void _JsonSetSequence(object set)
        {
            Sequence = (List < LadderOperator < TIndex >>) set;
        }

        /// <summary>
        /// Returns sign coefficient of ladder operator sequence.
        /// </summary>
        /// <returns>Sign of ladder operator sequence.</returns>
        public int _JsonGetCoefficient() => Coefficient;

        /// <summary>
        /// Sets sign coefficient of ladder operator sequence.
        /// </summary>
        /// <returns>Sign of ladder operator sequence.</returns>
        public void _JsonSetCoefficient(int set)
        {
            Coefficient = set;
        }
        /// <summary>
        /// Used only for JSON serialization.
        /// </summary>
        public void _JsonSetObject(object set)
        {
            var result = (ValueTuple<List<LadderOperator<TIndex>>, int>)set;
            Sequence = result.Item1;
            Coefficient = result.Item2;
        }
        #endregion

        #region Constructors
        /// <summary>
        /// Constructor for empty ladder operator sequence.
        /// </summary>
        public LadderSequence() : base() { }

        /// <summary>
        /// Constructor for an identitcal ladder operator sequence.
        /// </summary>
        public LadderSequence(LadderSequence<TIndex> setSequence)
        {
            Coefficient = setSequence.Coefficient;
            Sequence = setSequence.Sequence.ToList();
        }
        
        /// <summary>
        /// Construct instance from sequence of ladder operators.
        /// </summary>
        /// <param name="setSequence">Sequence of ladder operators.</param>
        public LadderSequence(IEnumerable<LadderOperator<TIndex>> setSequence, int setCoefficient = 1)
        {
            Sequence = setSequence.ToList();
            Coefficient = setCoefficient;
        }

        // This exists as a convenience function for creating fermion terms in samples.
        /// <summary>
        /// Implicit operator for creating a Ladder operator.
        /// </summary>
        /// <param name="setOperator">Tuple where the first parameter
        /// is the raising or lowering index, and the second parameter
        /// is the position index of the ladder operator.</param>
        public static implicit operator LadderSequence<TIndex>((RaisingLowering, TIndex)[] setSequence)
        {
            return new LadderSequence<TIndex>(setSequence.Select(o => new LadderOperator<TIndex>(o)));
        }

        /// <summary>
        /// Construct a sequence of ladder operators from an even-length sequence of integers.
        /// </summary>
        /// <param name="indices">Even-length sequence of integers.</param>
        /// <returns>
        /// Sequence of ladder operators with an equal number of creation and annihilation terms
        /// that are normal-ordered.
        /// </returns>
        /// <example>
        /// <code>
        /// // The following two return the same ladder operator sequence.
        /// var seq = new[] { 1, 2, 3, 4 }.ToLadderSequence();
        /// var expected = new[] { (u, 1), (u, 2), (d, 3), (d, 4) }.ToLadderSequence();
        /// </code>
        /// </example>
        public static implicit operator LadderSequence<TIndex>(TIndex[] indices)
        {
            var length = indices.Count();
            if (length % 2 == 1)
            {
                throw new System.ArgumentException(
                    $"Number of terms provided is `{length}` and must be of even length."
                    );
            }
            (RaisingLowering, TIndex) GetLadderOperator(TIndex index, int position)
            {
                return (position < length / 2 ? RaisingLowering.u : RaisingLowering.d, index);
            }
            return indices.Select(GetLadderOperator).ToArray();
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
        public bool IsInNormalOrder() => Sequence.Count() == 0 ? true : Sequence.Select(o => (int)o.Type).IsInAscendingOrder();
        #endregion

        /// <summary>
        /// Creates a new ladder sequence with a different indexing scheme.
        /// </summary>
        /// <typeparam name="TNewIndex">Type of the new indexing scheme.</typeparam>
        /// <param name="indexFunction">Function for mapping the current scheme to the new scheme.</param>
        /// <returns>Ladder sequence with a new index type.</returns>
        public LadderSequence<TNewIndex> SelectIndex<TNewIndex>(Func<TIndex, TNewIndex> indexFunction)
        where TNewIndex : IEquatable<TNewIndex>
        {
            var newIndexing = Sequence
                .Select(o => new LadderOperator<TNewIndex>
                (o.Type, indexFunction((TIndex) o.Index)));
            return new LadderSequence<TNewIndex>(newIndexing, Coefficient);
        }

        /// <summary>
        /// Concatenates two Fermion terms.
        /// </summary>
        /// <param name="left">Left <see cref="LadderSequence"/> <c>x</c>.</param>
        /// <param name="right">Right <see cref="LadderSequence"/> <c>y</c>.</param>
        /// <returns>
        /// Returns new <see cref="LadderSequence"/> <c>xy</c> where coefficients and 
        /// LadderOperatorSequences are multipled together.
        /// </returns>
        // TODO: May decide to overload the * operator.
        public LadderSequence<TIndex> Multiply(
            LadderSequence<TIndex> left, 
            LadderSequence<TIndex> right) 
            => new LadderSequence<TIndex>(left.Sequence.Concat(right.Sequence), left.Coefficient * right.Coefficient);
        
        /// <summary>
        /// Counts the number of unique system indices across all <see cref="LadderOperator"/> terms
        /// in a <see cref="LadderSequence"/>
        /// </summary>
        /// <returns>Number of unique system indices.</returns>
        public int UniqueIndices() => Sequence.Select(o => o.Index).Distinct().Count();

        /// <summary>
        /// Returns a copy of the ladder sequence base class.
        /// </summary>
        /// <returns>Base class of this sequence of ladder operators.</returns>
        public LadderSequence<TIndex> ToLadderSequence() => new LadderSequence<TIndex>(Sequence, Coefficient);

        /// <summary>
        /// Returns list of indices of the ladder operator sequence.
        /// </summary>
        /// <returns>Sequence of integers. </returns>
        public IEnumerable<TIndex> ToIndices() => Sequence.Select(o => (TIndex) o.Index);

        /// <summary>
        /// Returns sequence of raising and lowering types of the ladder operator sequence.
        /// </summary>
        /// <returns>Sequence of raising an lowering operators.</returns>
        public IEnumerable<RaisingLowering> ToRaisingLowering() => Sequence.Select(o => o.Type);

        /// <summary>
        /// Returns a human-readable description of this object.
        /// </summary>
        public override string ToString() => $"{Coefficient} * {string.Join(" ", Sequence)}";

        #region Equality Testing
        public override bool Equals(object obj) => (obj is LadderSequence<TIndex> x) ? Equals(x) : false;

        public bool Equals(LadderSequence<TIndex> x)
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

        public static bool operator ==(LadderSequence<TIndex> x, LadderSequence<TIndex> y)
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
        
        public static bool operator !=(LadderSequence<TIndex> x, LadderSequence<TIndex> y) => !(x == y);

        #endregion
    }
}



