// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using Microsoft.Quantum.Simulation.Core;

using System;
using Microsoft.Quantum.Chemistry.Generic;
using Microsoft.Quantum.Chemistry.Fermion;

namespace Microsoft.Quantum.Chemistry.QSharpFormat
{
    using Microsoft.Quantum.Chemistry.JordanWigner;
    
    /// <summary>
    /// Methods for converting electronic structure problem to data for consumption by Q#.
    /// </summary>
    public static partial class Convert
    {
        /// <summary>
        /// Translate initial state to a format consumable by Q#.
        /// </summary>
        /// <param name="inputState">Initial state</param>
        /// <returns>Initial state in Q# format.</returns>
        public static JordanWignerEncodingData ToQSharpFormat(
            (Double, Int64, JWOptimizedHTerms) pauliHamiltonianQSharpFormat,
            (Int64, QArray<JordanWignerInputState>) wavefunctionQSharpFormat
            )
        {
            var energyOffset = pauliHamiltonianQSharpFormat.Item1;
            var nSpinOrbitals = pauliHamiltonianQSharpFormat.Item2;

            return new JordanWignerEncodingData(
                (nSpinOrbitals
                , pauliHamiltonianQSharpFormat.Item3
                , wavefunctionQSharpFormat
                , energyOffset));
        }
    }

}