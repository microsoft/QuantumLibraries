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
    /// This JsonConverter encodes of a LadderOperator as a Tuple instead of as an object.
    /// </summary>
    public class LadderOperatorJsonConverter : JsonConverter
    {
        public override bool CanConvert(Type objectType)
        {
            return objectType.IsGenericType && (objectType.GetGenericTypeDefinition() == typeof(LadderOperator<>));
        }

        /// <summary>
        /// Writers the LadderOperator as a (Type, Index) tuple.
        /// </summary>
        public override void WriteJson(JsonWriter writer, object value, JsonSerializer serializer)
        {
            var thing = (ILadderOperator)value;
            var item = (thing.GetRaisingLowering(), thing.ObjectGetIndex());
            serializer.Serialize(writer, item);
        }

        /// <summary>
        /// Reads the LadderOperator from a (Type, Index) tuple.
        /// </summary>
        public override object ReadJson(JsonReader reader, Type objectType, object existingValue, JsonSerializer serializer)
        {
            var indexType = objectType.GetGenericArguments()[0];

            var tupleType = typeof(ValueTuple<,>).MakeGenericType(typeof(RaisingLowering), indexType);
            var item = serializer.Deserialize(reader, tupleType);

            var ladderOperatorType = typeof(LadderOperator<>).MakeGenericType(indexType);
            var result = (ILadderOperator)Activator.CreateInstance(ladderOperatorType);

            result.SetObject(item);

            return result;
        }
    }

    /// <summary>
    /// Ladder operators implement this interface. This interface is used to
    /// enable generic Json serialization.
    /// </summary>
    /// <typeparam name="TIndex">Ladder operator index type.</typeparam>
    public interface ILadderOperator
    {
        void SetObject(object set);
        object ObjectGetIndex();
        void SetIndex(object set);
        RaisingLowering GetRaisingLowering();
        void SetRaisingLowering(object set);

    }
}