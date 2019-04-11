// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using System;
using System.Linq;
using System.Collections.Generic;
using Newtonsoft.Json;

namespace Microsoft.Quantum.Chemistry.LadderOperators
{

    /// <summary>
    /// Class representing a sequence of raising and lowering operators.
    /// </summary>
    [JsonConverter(typeof(LadderSequence.JsonConverter))]
    public class LadderSequence : IEquatable<LadderSequence>
    {

        /// <summary>
        /// Sequence of ladder operators.
        /// </summary>
        public List<LadderOperator> Sequence { get; set; } = new List<LadderOperator>();

        /// <summary>
        /// sign (-1,+1) coefficient of ladder operators.
        /// </summary>
        public int Coefficient { get; set; } = 1;

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
            if (ladderOperators != null)
            {
                // All constructions are pass by value. 
                // ToList() copies values if the underlying object is a value type.
                Sequence = ladderOperators.Sequence.ToList();
                Coefficient = ladderOperators.Coefficient;
            }
        }

        /// <summary>
        /// Construct a sequence of ladder operators from sequence of tuples each
        /// specifying whether it is a raising or lowering term, and its index.
        /// </summary>
        /// <param name="setSequence">Sequence of ladder operators.</param>
        /// <param name="setSign">Set the sign coefficient of the sequence.</param>
        /// <returns>
        /// Sequence of ladder operators.
        /// </returns>
        /// <example>
        /// // Construct a sequence a ladder operators 1^ 2^ 3 4
        /// var tmp = new[] { (u, 1), (u, 2), (d, 3), (d, 4) }.ToLadderSequence();
        /// </example>
        public static implicit operator LadderSequence((RaisingLowering, int)[] setSequence)
        {
            return new LadderSequence(setSequence.Select(o => new LadderOperator(o)));
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
        public static implicit operator LadderSequence(int[] indices)
        {
            var length = indices.Count();
            if (length % 2 == 1)
            {
                throw new System.ArgumentException(
                    $"Number of terms provided is `{length}` and must be of even length."
                    );
            }
            Func<int, int, (RaisingLowering, int)> GetLadderOperator = (index, position)
                => (position < length / 2 ? RaisingLowering.u : RaisingLowering.d, index);
            return indices.Select((o, idx) => GetLadderOperator(o, idx)).ToLadderSequence();
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
        public bool IsInNormalOrder() => Sequence.Count() == 0 ? true : Sequence.Select(o => (int)o.Type).IsInAscendingOrder();
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
        // TODO: May decide to overload the * operator.
        public LadderSequence Multiply(LadderSequence left, LadderSequence right) => new LadderSequence(left.Sequence.Concat(right.Sequence), left.Coefficient * right.Coefficient);

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
        public int UniqueIndices() => Sequence.Select(o => o.Index).Distinct().Count();

        /// <summary>
        /// Returns a copy of the ladder sequence base class.
        /// </summary>
        /// <returns>Base class of this sequence of ladder operators.</returns>
        public LadderSequence GetLadderSequence() => new LadderSequence(Sequence); 

        /// <summary>
        /// Returns list of indices of the ladder operator sequence.
        /// </summary>
        /// <returns>Sequence of integers. </returns>
        public IEnumerable<int> Indices() => Sequence.Select(o => o.Index);

        /// <summary>
        /// Returns sequence of raising and lowering types of the ladder operator sequence.
        /// </summary>
        /// <returns>Sequence of raising an lowering operators.</returns>
        public IEnumerable<RaisingLowering> SequenceRaisingLowering() =>  Sequence.Select(o => o.Type);

        /// <summary>
        /// Returns a human-readable description this object.
        /// </summary>
        public override string ToString() => $"{Coefficient} * {string.Join(" ",Sequence)}";


        #region Equality Testing
        public override bool Equals(object obj) => (obj is LadderSequence x) ? Equals(x) : false;

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

        public static bool operator !=(LadderSequence x, LadderSequence y) => !(x == y);

        #endregion


        #region JsonConverter
        /// <summary>
        /// This JsonConverter encodes of a LadderSequence as a Tuple instead of as an object.
        /// </summary>
        public class JsonConverter : JsonConverter<LadderSequence>
        {
            /// <summary>
            /// Writers the LadderSequence as a (Sequence, Coefficient) tuple.
            /// </summary>
            public override void WriteJson(JsonWriter writer, LadderSequence value, JsonSerializer serializer)
            {
                if (value == null)
                {
                    serializer.Serialize(writer, null);
                }
                else
                {
                    var item = (value.Sequence, value.Coefficient);
                    serializer.Serialize(writer, item);
                }
            }

            /// <summary>
            /// Reads the LadderSequence from aa (Sequence, Coefficient) tuple.
            /// </summary>
            public override LadderSequence ReadJson(JsonReader reader, Type objectType, LadderSequence existingValue, bool hasExistingValue, JsonSerializer serializer)
            {
                if (reader.TokenType == JsonToken.Null)
                    return null;

                var item = serializer.Deserialize<ValueTuple<List<LadderOperator>, int>>(reader);

                var result = Activator.CreateInstance(objectType) as LadderSequence;
                result.Sequence = item.Item1;
                result.Coefficient = item.Item2;

                return result;
            }
        }
        #endregion
    }
}



