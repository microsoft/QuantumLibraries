// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using System;
using Newtonsoft.Json;

namespace Microsoft.Quantum.Chemistry.LadderOperators
{


    /// <summary>
    /// Data strcture for raising and lowering operators.
    /// </summary>
    [JsonConverter(typeof(LadderOperator.JsonConverter))]
    public struct LadderOperator : IEquatable<LadderOperator>
    {
        /// <summary>
        /// LadderType specifying raising or lowering operator.
        /// </summary>
        public RaisingLowering Type;

        /// <summary>
        /// System index operator acts on.
        /// </summary>
        public int Index;
        

        /// <summary>
        /// Constructor for ladder operator.
        /// </summary>
        /// <param name="setType">Set raising or lowering operator.</param>
        /// <param name="setIndex">Set system index.</param>
        public LadderOperator(RaisingLowering setType, int setIndex)
        {
            Type = setType;
            Index = setIndex;
        }

        /// <summary>
        /// Constructor for ladder operator.
        /// </summary>
        /// <param name="set">Tuple where first item sets the operator type,
        /// and the second item indexes the system.</param>
        public LadderOperator((RaisingLowering, int) set) : this(set.Item1, set.Item2) { }

        #region Equality Testing
        public override bool Equals(object obj) => (obj is LadderOperator x) ? Equals(x) : false;

        public bool Equals(LadderOperator x) => (Type == x.Type) && (Index == x.Index);

        public override int GetHashCode() => Index * Enum.GetNames(typeof(RaisingLowering)).Length + (int) Type;

        public static bool operator ==(LadderOperator x, LadderOperator y) => x.Equals(y);

        public static bool operator !=(LadderOperator x, LadderOperator y) => !(x.Equals(y));
        #endregion

        /// <summary>
        /// Returns a human-readable description this object.
        /// </summary>
        public override string ToString() {
            string op = Type.ToString();
            return $"{Index}{op}";
        }


        /// <summary>
        /// This JsonConverter encodes of a LadderOperator as a Tuple instead of as an object.
        /// </summary>
        public class JsonConverter : JsonConverter<LadderOperator>
        {
            /// <summary>
            /// Writers the LadderOperator as a (Type, Index) tuple.
            /// </summary>
            public override void WriteJson(JsonWriter writer, LadderOperator value, JsonSerializer serializer)
            {
                var item = (value.Type, value.Index);
                serializer.Serialize(writer, item);
            }

            /// <summary>
            /// Reads the LadderOperator from a (Type, Index) tuple.
            /// </summary>
            public override LadderOperator ReadJson(JsonReader reader, Type objectType, LadderOperator existingValue, bool hasExistingValue, JsonSerializer serializer)
            {
                var item = serializer.Deserialize<ValueTuple<RaisingLowering, int>>(reader);

                var result = (LadderOperator)Activator.CreateInstance(objectType);
                result.Type = item.Item1;
                result.Index = item.Item2;

                return result;
            }
        }
    }
    

}



