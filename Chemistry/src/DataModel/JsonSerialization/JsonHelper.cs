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
  
    internal static class Helper
    {
        internal static Type GetBasestType(Type t)
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

}