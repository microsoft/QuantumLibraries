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
    public static class Wavefunction
    {
        [Fact]
        static void MakeHartreeFock()
        {
            // We initialize a fermion Hamiltonian.
            var fermionHamiltonian = new FermionHamiltonian();

            // Create a Hartree-Fock state from the Hamiltonian 
            // with, say, `4` occupied spin orbitals.
            var wavefunction = fermionHamiltonian.CreateHartreeFockState(nElectrons: 4);
        }
    }

}
