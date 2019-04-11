using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;
using System.Text;
using Newtonsoft.Json;

namespace Microsoft.Quantum.Chemistry.LadderOperators
{
    /// <summary>
    /// This JsonConverter encodes of a LadderSequence as a Tuple instead of as an object.
    /// </summary>
    public class LadderSequenceJsonConverter : JsonConverter<LadderSequence>
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

            var item = serializer.Deserialize<ValueTuple<List<LadderOperator>,int>>(reader);

            var result = Activator.CreateInstance(objectType) as LadderSequence;
            result.Sequence = item.Item1;
            result.Coefficient = item.Item2;

            return result;
        }
    }

    /// <summary>
    /// This JsonConverter encodes of a LadderOperator as a Tuple instead of as an object.
    /// </summary>
    public class LadderOperatorJsonConverter : JsonConverter<LadderOperator>
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