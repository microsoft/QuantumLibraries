using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;
using System.Text;
using Newtonsoft.Json;

namespace Microsoft.Quantum.Chemistry.Generic
{
    /// <summary>
    /// This JsonConverter allows to correctly serialized HamiltonianTerms.
    /// This terms are in general problematic because their keys are not strings, 
    /// but HamiltonianTerms, which json.net doesn't like by default. 
    /// This converts the Dictionaries to List of Tuples, in which the first
    /// item of the tuple is the key and the second the value.
    /// </summary>
    public class HamiltonianTermsJsonConverter : JsonConverter
    {
        /// <summary>
        /// Returns true only if the Type is HamitonianTerm or HamiltonianTerms
        /// </summary>
        public override bool CanConvert(Type objectType) =>
            objectType.IsGenericType
            && (objectType.GetGenericTypeDefinition() == typeof(Hamiltonian<,,>.HamiltonianTerm)
                || objectType.GetGenericTypeDefinition() == typeof(Hamiltonian<,,>.HamiltonianTerms));

        /// <summary>
        /// Writers the HamiltonianTerms as a list of (Key, Value) tuples.
        /// </summary>
        public override void WriteJson(JsonWriter writer, object value, JsonSerializer serializer)
        {
            var dictionary = (IDictionary)value;

            writer.WriteStartArray();
            foreach (var key in dictionary.Keys)
            {
                var item = (key, dictionary[key]);
                serializer.Serialize(writer, item);
            }
            writer.WriteEndArray();
        }

        /// <summary>
        /// Reads the HamiltonianTerms from a list of (Key, Value) tuples.
        /// </summary>
        public override object ReadJson(JsonReader reader, Type objectType, object existingValue, JsonSerializer serializer)
        {
            Debug.Assert(CanConvert(objectType));

            var keyType = objectType.BaseType.GetGenericArguments()[0];
            var valueType = objectType.BaseType.GetGenericArguments()[1];
            var tupleType = typeof(ValueTuple<,>).MakeGenericType(keyType, valueType);
            
            var result = (IDictionary)Activator.CreateInstance(objectType);

            if (reader.TokenType == JsonToken.Null)
                return null;

            while (reader.Read())
            {
                if (reader.TokenType == JsonToken.EndArray)
                {
                    return result;
                }

                if (reader.TokenType == JsonToken.StartObject)
                {
                    dynamic item = serializer.Deserialize(reader, tupleType);
                    result.Add(item.Item1, item.Item2);
                }
            }

            return result;
        }
    }
}