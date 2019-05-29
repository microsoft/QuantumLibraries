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
        /// Combine Q# Hamiltonian and wavefunciton data to a format consumed
        /// by the Q# chemistry libraries.
        /// </summary>
        /// <param name="pauliHamiltonianQSharpFormat">Hamiltonian data in Q# format.</param>
        /// <param name="wavefunctionQSharpFormat">Wavefunction data in Q# format.</param>
        /// <returns>Combined Hamiltonian and wave function data in Q# format.</returns>
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

    /// <summary>
    /// Methods for converting electronic structure problem to data for consumption by Q#.
    /// </summary>
    public static partial class Extension
    {
        /// <summary>
        /// Combine Q# Hamiltonian with empty wavefunciton data to a format consumed
        /// by the Q# chemistry libraries.
        /// </summary>
        /// <param name="pauliHamiltonianQSharpFormat">Hamiltonian data in Q# format.</param>
        /// <returns>Combined Hamiltonian and wave function data in Q# format.</returns>
        public static JordanWignerEncodingData Pad(
            this (Double, Int64, JWOptimizedHTerms) pauliHamiltonianQSharpFormat
            )
        {
            (Int64, QArray<JordanWignerInputState>) wavefunctionQSharpFormat = (0, new QArray<JordanWignerInputState>());

            return Convert.ToQSharpFormat(pauliHamiltonianQSharpFormat, wavefunctionQSharpFormat);
        }

        /// <summary>
        /// Combine Q# wavefunction data with empty Hamiltonian data to a format consumed
        /// by the Q# chemistry libraries.
        /// </summary>
        /// <param name="wavefunctionQSharpFormat">Wavefunction data in Q# format.</param>
        /// <returns>Combined Hamiltonian and wave function data in Q# format.</returns>
        public static JordanWignerEncodingData Pad(
            this (Int64, QArray<JordanWignerInputState>) wavefunctionQSharpFormat
            )
        {

            (Double, Int64, JWOptimizedHTerms) pauliHamiltonianQSharpFormat = (0.0, 0, new JWOptimizedHTerms());
            return Convert.ToQSharpFormat(pauliHamiltonianQSharpFormat, wavefunctionQSharpFormat);
        }
    }

}