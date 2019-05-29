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

namespace Microsoft.Quantum.Chemistry.Tests.Samples.Hydrogen
{
    public static class MolecularHydrogen
    {
        [Fact]
        static void MolecularHydrogenTest()
        {
            //////////////////////////////////////////////////////////////////////////
            // Introduction //////////////////////////////////////////////////////////
            //////////////////////////////////////////////////////////////////////////

            // In this example, we will create a spin-orbital representation of the molecular
            // Hydrogen Hamiltonian `H`, given ovelap coefficients for its one- and 
            // two - electron integrals.

            // We when perform quantum phase estimation to obtain an estimate of
            // the molecular Hydrogen ground state energy.

            #region Building the Hydrogen Hamiltonian through orbital integrals

            // One of the simplest representations of Hydrogen uses only two 
            // molecular orbitals indexed by `0` and `1`.
            var nOrbitals = 2;

            // This representation also has two occupied spin-orbitals.
            var nElectrons = 2;

            // The Coulomb repulsion energy between nuclei is
            var energyOffset = 0.713776188;

            // One-electron integrals are listed below
            // <0|H|0> = -1.252477495
            // <1|H|1> = -0.475934275

            // Two-electron integrals are listed below
            // <00|H|00> = 0.674493166
            // <01|H|01> = 0.181287518
            // <01|H|10> = 0.663472101
            // <11|H|11> = 0.697398010
            // Note that orbitals are assumed to be real. Moreover, indistinguishability
            // of electrons means that the following integrals are equal.
            //   <PQ|H|RS> = <PR|H|QS> = <SQ|H|RP> = <SR|H|QP>
            // = <QP|H|SR> = <RP|H|SQ> = <QS|H|PR> = <RS|H|PQ>
            // Thus it sufficies to specify just any one of these terms from each symmetry
            // group.

            // These orbital integrals are represented using the OrbitalIntegral
            // data structure.
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

            // These orbital integral terms are automatically expanded into
            // spin-orbitals. We may print the Hamiltonian to see verify what it contains.
            Console.WriteLine("----- Print Hamiltonian");
            Console.Write(fermionHamiltonian);
            Console.WriteLine("----- End Print Hamiltonian \n");

            // We also need to create an input quantum state to this Hamiltonian.
            // Let us use the Hartree-Fock state.
            var fermionWavefunction = fermionHamiltonian.CreateHartreeFockState(nElectrons);
            #endregion

            #region Jordan-Wigner representation 
            // The Jordan-Wigner encoding converts the Fermion Hamiltonian, 
            // expressed in terms of Fermionic operators, to a qubit Hamiltonian,
            // expressed in terms of Pauli matrices. This is an essential step
            // for simulating our constructed Hamiltonians on a qubit quantum
            // computer.
            Console.WriteLine("----- Creating Jordan-Wigner encoding");
            var jordanWignerEncoding = fermionHamiltonian.ToPauliHamiltonian(Pauli.QubitEncoding.JordanWigner);
            Console.WriteLine("----- End Creating Jordan-Wigner encoding \n");
            #endregion

            #region Performing the simulation 
            // We are now ready to run a quantum simulation of molecular Hydrogen.
            // We will use this to obtain an estimate of its ground state energy.


                // This Jordan-Wigner data structure also contains a representation 
                // of the Hamiltonian and wavefunction made for consumption by the Q# algorithms.
                var qSharpHamiltonianData = jordanWignerEncoding.ToQSharpFormat();
                var qSharpWavefunctionData = fermionWavefunction.ToQSharpFormat();
                var qSharpData = QSharpFormat.Convert.ToQSharpFormat(qSharpHamiltonianData, qSharpWavefunctionData);

                

            Console.WriteLine("Press Enter to continue...");
            if (System.Diagnostics.Debugger.IsAttached)
            {
                Console.ReadLine();
            }
            #endregion

        }
    }
}
