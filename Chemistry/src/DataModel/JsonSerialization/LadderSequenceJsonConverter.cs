using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;
using System.Text;
using Newtonsoft.Json;
using Microsoft.Quantum.Chemistry.OrbitalIntegrals;

using Microsoft.Quantum.Chemistry.LadderOperators;

// Main idea: We do not want to implement a general technique to serialize everything.
// We will create wrapper classes around the specific objects we want to serialize.


// Serialize Ladder operators: Int and SpinOrbital type.
// Serialize Hamiltonian: Fermion type.

   

namespace Microsoft.Quantum.Chemistry.Json
{
  
    /// <summary>
    /// This JsonConverter encodes of a LadderSequence as a Tuple instead of as an object.
    /// </summary>
    public class LadderSequenceJsonConverter : JsonConverter
    {
        public override bool CanConvert(Type objectType)
        {
            return objectType.IsGenericType && (objectType.GetGenericTypeDefinition() == typeof(LadderSequence<>));
        }
        /// <summary>
        /// Writers the LadderSequence as a (Sequence, Coefficient) tuple.
        /// </summary>
        public override void WriteJson(JsonWriter writer, object value, JsonSerializer serializer)
        {
            if (value == null)
            {
                serializer.Serialize(writer, null);
            }
            else
            {
                var thing = (ILadderSequence)value;

                var item = (thing.GetSequence(), thing.GetCoefficient());
                serializer.Serialize(writer, item);
            }
        }

        /// <summary>
        /// Reads the LadderSequence from aa (Sequence, Coefficient) tuple.
        /// </summary>
        public override object ReadJson(JsonReader reader, Type objectType, object existingValue, JsonSerializer serializer)
        {
            if (reader.TokenType == JsonToken.Null)
            { 
                return null;
            }
            
            var indexType = GetBasestType(objectType).GetGenericArguments()[0];
            var tupleType = typeof(ValueTuple<,>).MakeGenericType( 
                typeof(List<>).MakeGenericType(
                    typeof(LadderOperator<>)
                    .MakeGenericType(indexType)), indexType);

            var item = serializer.Deserialize(reader, tupleType);

            var result = (ILadderSequence)Activator.CreateInstance(objectType);
            result.SetObject(item);

            return result;
        }

        internal Type GetBasestType(Type t)
        {
            var type = t;
            if (t.BaseType == typeof(object))
            {
                return t;
            }
            else
            {
                return GetBasestType(t.BaseType);
            }
        }
    }
    
    /// <summary>
    /// Ladder sequences implement this interface. This interface is used to
    /// enable generic Json serialization.
    /// </summary>
    /// <typeparam name="TIndex">Ladder sequence index type.</typeparam>
    public interface ILadderSequence
    {
        void SetObject(object set);
        object GetSequence();
        void SetSequence(object set);
        int GetCoefficient();
        void SetCoefficient(int set);

    }

}