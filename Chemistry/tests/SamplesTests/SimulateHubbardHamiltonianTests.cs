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
using Microsoft.Quantum.Chemistry;
using Microsoft.Quantum.Chemistry.OrbitalIntegrals;
using Microsoft.Quantum.Chemistry.Fermion;
using Microsoft.Quantum.Chemistry.QSharpFormat;

// To count gates, we'll use the trace simulator provided with
// the Quantum Development Kit.
using Microsoft.Quantum.Simulation.Simulators.QCTraceSimulators;

// The System namespace provides a number of useful built-in
// types and methods that we'll use throughout this sample.
using System;

// The System.Diagnostics namespace provides us with the
// Stopwatch class, which is quite useful for measuring
// how long each gate counting run takes.
using System.Diagnostics;

// The System.Collections.Generic library provides many different
// utilities for working with collections such as lists and dictionaries.
using System.Collections.Generic;

// We use the logging library provided with .NET Core to handle output
// in a robust way that makes it easy to turn on and off different messages.
using Microsoft.Extensions.Logging;

// We use this for convnience functions for manipulation arrays.
using System.Linq;

//
using Xunit;
#endregion

namespace Microsoft.Quantum.Chemistry.Tests.Samples.Hubbard
{
    public static class SimulateHubbardHamiltonian
    {
        [Fact]
        static void SimulateHubbardHamiltonianTest()
        {
            //////////////////////////////////////////////////////////////////////////
            // Introduction //////////////////////////////////////////////////////////
            //////////////////////////////////////////////////////////////////////////

            // In this example, we will estimate the ground state energy of 
            // 1D Hubbard Hamiltonian using the quantum chemistry library. 

            // The 1D Hubbard model has `n` sites. Let `i` be the site index, 
            // `s` = 1,0 be the spin index, where 0 is up and 1 is down, `t` be the 
            // hopping coefficient, `u` the repulsion coefficient, and aᵢₛ the fermionic 
            // annihilation operator on the fermion indexed by `(i,s)`. The Hamiltonian 
            // of this model is
            //
            //     H ≔ - t Σᵢ (a†ᵢₛ aᵢ₊₁ₛ + a†ᵢ₊₁ₛ aᵢₛ) + u Σᵢ a†ᵢ₀ a†ᵢ₁ aᵢ₁ aᵢ₀
            //
            // Note that we use closed boundary conditions.

            #region Building the Hubbard Hamiltonian through orbital integrals

            var t = 0.2; // hopping coefficient
            var u = 1.0; // repulsion coefficient
            var nSites = 6; // number of sites;
            // Construct Hubbard Hamiltonian
            var hubbardOrbitalIntegralHamiltonian = new OrbitalIntegralHamiltonian();

            foreach (var i in Enumerable.Range(0, nSites))
            {
                hubbardOrbitalIntegralHamiltonian.Add(new OrbitalIntegral(new[] { i, (i + 1) % nSites }, -t));
                hubbardOrbitalIntegralHamiltonian.Add(new OrbitalIntegral(new[] { i, i, i, i }, u));
            }

            // Create fermion representation of Hamiltonian
            // In this case, we use the spin-orbital to integer
            // indexing convention `x = orbitalIdx + spin * nSites`; as it 
            // minimizes the length of Jordan–Wigner strings
            var hubbardFermionHamiltonian = hubbardOrbitalIntegralHamiltonian.ToFermionHamiltonian(IndexConvention.HalfUp);

            #endregion


            #region Estimating energies by simulating quantum phase estimation
            // Create Jordan–Wigner representation of Hamiltonian
            var jordanWignerEncoding = hubbardFermionHamiltonian.ToPauliHamiltonian();

            // Create data structure to pass to QSharp.
            var qSharpData = jordanWignerEncoding.ToQSharpFormat().Pad();

            Console.WriteLine($"Estimate Hubbard Hamiltonian energy:");

            #endregion
        }
    }
}
