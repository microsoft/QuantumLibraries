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
  
    internal static class TypeExtensions
    {
        /// <summary>
        ///       Searches base types of a given type to find the type that immediately derives from
        ///       <see href="System.Object" />.
        /// </summary>
        internal static Type GetBasestType(this Type t) =>
            (t.BaseType == typeof(object) || t == typeof(object)) ? t : GetBasestType(t.BaseType);
    }

}