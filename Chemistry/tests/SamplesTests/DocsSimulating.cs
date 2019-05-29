// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

// This test ensures that any chemistry library syntax changes
// that affect samples are detected.

#region Using Statements
// We will need several different libraries in this sample.
// Here, we expose these libraries to our program using the
// C# "using" statement, similar to the Q# "open" statement.

// We will use the data model implemented by the Quantum Development Kit Chemistry
// Libraries. This model defines what a fermionic Hamiltonian is, and how to
// represent Hamiltonians on disk.
using Microsoft.Quantum.Chemistry.OrbitalIntegrals;
using Microsoft.Quantum.Chemistry.Fermion;
using Microsoft.Quantum.Chemistry.Paulis;
using Microsoft.Quantum.Chemistry.QSharpFormat;

// The System namespace provides a number of useful built-in
// types and methods that we'll use throughout this sample.
using System;

// We use this for convnience functions for manipulation arrays.
using System.Linq;

//
using Xunit;
#endregion

namespace Microsoft.Quantum.Chemistry.Tests.Docs
{
    public class SimulatingHamiltonianDynamics
    {
        [Fact]
        public void MakeHamiltonian()
        {
            // We create a `FermionHamiltonian` object to store the terms.
            var hamiltonian = new OrbitalIntegralHamiltonian(
                new[] {
                new OrbitalIntegral(new[] { 0, 1, 2, 3 }, 0.123),
                new OrbitalIntegral(new[] { 0, 1 }, 0.456)
                    }
                )
                .ToFermionHamiltonian(IndexConvention.UpDown);

            // We convert this fermion Hamiltonian to a Jordan-Wigner representation.
            var jordanWignerEncoding = hamiltonian.ToPauliHamiltonian(QubitEncoding.JordanWigner);

            // We now convert this representation into a format consumable by Q#.
            var qSharpData = jordanWignerEncoding.ToQSharpFormat();

            Assert.True(hamiltonian.CountTerms() == 10);
            Assert.True(jordanWignerEncoding.CountTerms() == 6);
        }
    }
}

