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
    public static class InvokingChemistryLibrary
    {
        [Fact]
        static void MakeHamiltoniam()
        {
            // These orbital integrals are represented using the OrbitalIntegral
            // data structure.
            var energyOffset = 0.713776188;
            var nElectrons = 2;
            var orbitalIntegrals = new OrbitalIntegral[]
            {
                new OrbitalIntegral(new[] { 0,0 }, -1.252477495),
                new OrbitalIntegral(new[] { 1,1 }, -0.475934275),
                new OrbitalIntegral(new[] { 0,0,0,0 }, 0.674493166),
                new OrbitalIntegral(new[] { 0,1,0,1 }, 0.181287518),
                new OrbitalIntegral(new[] { 0,1,1,0 }, 0.663472101),
                new OrbitalIntegral(new[] { 1,1,1,1 }, 0.697398010),
                // Add the identity term
                new OrbitalIntegral(new int[] { }, energyOffset)
            };

            // We initialize a fermion Hamiltonian data structure and add terms to it
            var fermionHamiltonian = new OrbitalIntegralHamiltonian(orbitalIntegrals).ToFermionHamiltonian();

            // The Jordan-Wigner encoding converts the Fermion Hamiltonian, 
            // expressed in terms of Fermionic operators, to a qubit Hamiltonian,
            // expressed in terms of Pauli matrices. This is an essential step
            // for simulating our constructed Hamiltonians on a qubit quantum
            // computer.
            var jordanWignerEncoding = fermionHamiltonian.ToPauliHamiltonian(Pauli.QubitEncoding.JordanWigner);

            // We also need to create an input quantum state to this Hamiltonian.
            // Let us use the Hartree-Fock state.
            var fermionWavefunction = fermionHamiltonian.CreateHartreeFockState(nElectrons);

            // This Jordan-Wigner data structure also contains a representation 
            // of the Hamiltonian and wavefunction made for consumption by the Q# algorithms.
            var qSharpHamiltonianData = jordanWignerEncoding.ToQSharpFormat();
            var qSharpWavefunctionData = fermionWavefunction.ToQSharpFormat();
            var qSharpData = QSharpFormat.Convert.ToQSharpFormat(qSharpHamiltonianData, qSharpWavefunctionData);
        }

        [Fact]
        static void LoadFromBroombridgeFile()
        {
            // This is the name of the file we want to load
            var filename = @"n2_1_00Re_sto3g.nw.out.yaml";
            // This is the directory containing the file
            var root = @"IntegralData\YAML";

            // This deserializes Broombridge.
            var broombridge = Broombridge.Deserializers.DeserializeBroombridge($@"{root}\{filename}");

            // Note that the deserializer returns a list of `ProblemDescriptions` instances 
            // as the file might describe multiple Hamiltonians. In this example, there is 
            // only one Hamiltonian. So we use `.First()`, which selects the first element of the list.
            var problem = broombridge.ProblemDescriptions.First();

            // This extracts the `OrbitalIntegralHamiltonian` from Broombridge format,
            // then converts it to a fermion Hamiltonian, then to a Jordan-Wigner
            // representation.
            var orbitalIntegralHamiltonian = problem.OrbitalIntegralHamiltonian;
            var fermionHamiltonian = orbitalIntegralHamiltonian.ToFermionHamiltonian(IndexConvention.UpDown);
            var jordanWignerEncoding = fermionHamiltonian.ToPauliHamiltonian(Pauli.QubitEncoding.JordanWigner);

            // The desired initial state, assuming that a description of it is present in the
            // Broombridge schema.
            var state = "|E1>";
            var wavefunction = problem.Wavefunctions[state].ToIndexing(IndexConvention.UpDown);

            // This creates the qSharpData consumable by the Q# chemistry library algorithms.
            var qSharpHamiltonianData = jordanWignerEncoding.ToQSharpFormat();
            var qSharpWavefunctionData = wavefunction.ToQSharpFormat();
            var qSharpData = QSharpFormat.Convert.ToQSharpFormat(qSharpHamiltonianData, qSharpWavefunctionData);

        }


        [Fact]
        static void LoadFromLiquidFile()
        {
            // This is the name of the file we want to load
            var filename = @"fe2s2_sto3g.dat"; // This is Ferrodoxin.
            // This is the directory containing the file
            var root = @"IntegralData\Liquid";

            // Deserialize the LiQuiD format.
            var problem = LiQuiD.Deserialize($@"{root}\{filename}").First();

            // This extracts the `OrbitalIntegralHamiltonian` from problem
            // description format.
            var orbitalIntegralHamiltonian = problem.OrbitalIntegralHamiltonian;
        }

        [Fact]
        static void ResourceEstimate()
        {
            // Filename of Hamiltonian to be loaded.
            var filename = "...";

            // This deserializes Broombridge.
            var problem = Broombridge.Deserializers.DeserializeBroombridge(filename).ProblemDescriptions.First();

            // This is a data structure representing the Jordan-Wigner encoding 
            // of the Hamiltonian that we may pass to a Q# algorithm.
            var qSharpData = problem.ToQSharpFormat();
        }

        [Fact]
        static void AbbreviatedResourceEstimate()
        {
            // Filename of Hamiltonian to be loaded.
            var filename = "...";

            // This deserializes Broombridge.
            var problem = Broombridge.Deserializers.DeserializeBroombridge(filename).ProblemDescriptions.First();

            // This is a data structure representing the Jordan-Wigner encoding 
            // of the Hamiltonian that we may pass to a Q# algorithm.
            var qSharpData = problem.ToQSharpFormat();
        }
    }

}
