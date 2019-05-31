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
    /// This JsonConverter encodes of a LadderOperator as a System.ValueTuple instead of as an object.
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
            var item = (thing._JsonGetRaisingLowering(), thing._JsonObjectGetIndex());
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

            result._JsonSetObject(item);

            return result;
        }
    }

    /// <summary>
    /// Ladder operators implement this interface. This interface is used to
    /// enable generic Json serialization.
    /// </summary>
    internal interface ILadderOperator
    {
        /// <summary>
        /// Sets the ladder operator parameters to be an instance 
        /// represented by this `object`.
        /// </summary>
        /// <param name="set">Ladder operator parameters settings.</param>
        void _JsonSetObject(object set);
        /// <summary>
        /// Returns the index of this ladder operator as an `object`.
        /// </summary>
        object _JsonObjectGetIndex();
        /// <summary>
        /// Sets the ladder operator index to be an instance 
        /// represented by this `object`.
        /// </summary>
        /// <param name="set">Ladder operator index settings.</param>
        void _JsonSetIndex(object set);
        /// <summary>
        /// Returns the type of this ladder operator.
        /// </summary>
        RaisingLowering _JsonGetRaisingLowering();
        /// <summary>
        /// Sets the ladder operator type to be an instance 
        /// represented by this `object`.
        /// </summary>
        /// <param name="set">Ladder operator type settings.</param>
        void _JsonSetRaisingLowering(object set);

    }
}