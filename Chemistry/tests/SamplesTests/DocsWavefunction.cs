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
using Microsoft.Quantum.Chemistry.LadderOperators;

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
        static void MakeSingleReferenceState()
        {
            // Create a list of indices of the creation operators
            var indices = new[] { 1, 2, 6 };

            // Convert the list of indices to a `FermionWavefunction` object.
            var wavefunction = new FermionWavefunction<int>(indices);
        }

        [Fact]
        static void MakeSingleReferenceStateSpinOrbital()
        {
            // Create a list of indices of the creation operators
            var indices = new[] { (0, Spin.d), (1,Spin.u), (3, Spin.u) };

            // Convert the list of indices to a `FermionWavefunction` object.
            var wavefunctionSpinOrbital = new FermionWavefunction<SpinOrbital>(indices.ToSpinOrbitals());

            // Convert a wavefunction indexed by spin orbitals to
            // one indexed by integers
            var wavefunctionInt = wavefunctionSpinOrbital.ToIndexing(IndexConvention.UpDown);
        }

        [Fact]
        static void MakeHartreeFock()
        {
            // We initialize a fermion Hamiltonian.
            var fermionHamiltonian = new FermionHamiltonian();

            // Create a Hartree-Fock state from the Hamiltonian 
            // with, say, `4` occupied spin orbitals.
            var wavefunction = fermionHamiltonian.CreateHartreeFockState(nElectrons: 4);
        }

        [Fact]
        static void MakeMCFState()
        {
            // Create a list of tuples where the first item of each 
            // tuple are indices to the creation operators acting on the
            // vacuum state, and the second item is the coefficient
            // of that basis state.
            var superposition = new[] {
                (new[] {1, 2, 6}, 0.1),
                (new[] {2, 1, 5}, -0.2) };

            // Create a fermion wavefunction object that represents the superposition.
            var wavefunction = new FermionWavefunction<int>(superposition);
        }

        [Fact]
        static void MakeUCCState()
        {
            // Create a list of indices of the creation operators
            // for the single-reference state
            var reference = new[] { 1, 2 };

            // Create a list describing the cluster operator
            // The first half of each list of integers will be
            // associated with the creation operators, and
            // the second half with the annihilation operators.
            var clusterOperator = new[]
            {
                (new [] {0, 1}, 0.123),
                (new [] {0, 3, 1, 2}, 0.456),
                (new [] {3, 2, 1, 0}, 0.789)
            };
            
            // Create a fermion wavefunction object that represents the 
            // unitary coupled-cluster wavefunction. It is assumed implicity
            // that the exponent of the unitary coupled-cluster operator
            // is the cluster operator minus its Hermitian conjugate.
            var wavefunction = new FermionWavefunction<int>(reference, clusterOperator);
        }

        [Fact]
        static void MakeUCCSpinOrbitalState()
        {
            // Create a list of indices of the creation operators
            // for the single-reference state
            var reference = new[] { (1, Spin.u), (2, Spin.d) }.ToSpinOrbitals();

            // Create a list describing the cluster operator
            // The first half of each list of integers will be
            // associated with the creation operators, and
            // the second half with the annihilation operators.
            var clusterOperator = new[]
            {
                (new [] {(0, Spin.u), (1, Spin.u)}, 0.123),
                (new [] {(0, Spin.u), (3, Spin.d), (1, Spin.u), (2, Spin.d)}, 0.456),
                (new [] {(3, Spin.u), (2, Spin.u), (1, Spin.u), (0, Spin.u)}, 0.789)
            }.Select(o => (o.Item1.ToSpinOrbitals(), o.Item2));

            // Create a fermion wavefunction object that represents the 
            // unitary coupled-cluster wavefunction. It is assumed implicity
            // that the exponent of the unitary coupled-cluster operator
            // is the cluster operator minus its Hermitian conjugate.
            var wavefunctionSpinOrbital = new FermionWavefunction<SpinOrbital>(reference, clusterOperator);

            // Convert the wavefunction indexed by spin-orbitals to one indexed
            // by integers
            var wavefunctionInteger = wavefunctionSpinOrbital.ToIndexing(IndexConvention.UpDown);
        }

        [Fact]
        static void MakeUCCAllExcitations()
        {
            // Create a list of indices of the creation operators
            // for the single-reference state
            var reference = new[] { (1, Spin.u), (2, Spin.d) }.ToSpinOrbitals();

            // Generate all spin-conversing excitations from spin-orbitals 
            // occupied by the reference state to virtual orbitals.
            var generatedExcitations = reference.CreateAllUCCSDSingletExcitations(nOrbitals: 3).Excitations;

            // This is the list of expected spin-consvering excitations
            var expectedExcitations = new[]
            {
                new []{ (0, Spin.u), (1,Spin.u)},
                new []{ (2, Spin.u), (1,Spin.u)},
                new []{ (0, Spin.d), (2,Spin.d)},
                new []{ (1, Spin.d), (2,Spin.d)},
                new []{ (0, Spin.u), (0, Spin.d), (2, Spin.d), (1,Spin.u)},
                new []{ (0, Spin.u), (1, Spin.d), (2, Spin.d), (1,Spin.u)},
                new []{ (0, Spin.d), (2, Spin.u), (2, Spin.d), (1,Spin.u)},
                new []{ (1, Spin.d), (2, Spin.u), (2, Spin.d), (1,Spin.u)}
            }.Select(o => new IndexOrderedSequence<SpinOrbital>(o.ToLadderSequence()));

            // The following two assertions are true, and verify that the generated 
            // excitations exactly match the expected excitations.
            var bool0 = generatedExcitations.Keys.All(expectedExcitations.Contains);
            var bool1 = generatedExcitations.Count() == expectedExcitations.Count();
        }
    }

}
