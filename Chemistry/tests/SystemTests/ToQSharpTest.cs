// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using Xunit;
using Microsoft.Quantum.Chemistry;
using Microsoft.Quantum.Simulation.Core;

using System.Text.RegularExpressions;
using System;
using System.Linq;
using System.Collections.Generic;

using Microsoft.Quantum.Chemistry.OrbitalIntegrals;
using Microsoft.Quantum.Chemistry.Fermion;
using Microsoft.Quantum.Chemistry.LadderOperators;
using Microsoft.Quantum.Chemistry.Pauli;
using Microsoft.Quantum.Chemistry.QSharpFormat;

namespace Microsoft.Quantum.Chemistry.Tests
{
    using Microsoft.Quantum.Chemistry.QSharpFormat;


    using FermionTerm = FermionTerm;
    using SpinOrbital = SpinOrbital;

    public class JordanWignerOptimizedEvolutionSetTests
    {
        public FermionHamiltonian GenerateTestHamiltonian()
        {
            var hamiltonian = new FermionHamiltonian();
            

            var nOrbitals = 6;
            var nElectrons = 2;

            hamiltonian.Add(new HermitianFermionTerm( new [] { 0, 0 }.ToLadderSequence()), 1.0);
            hamiltonian.Add(new HermitianFermionTerm( new [] { 1, 1 }.ToLadderSequence()), 1.0);
            hamiltonian.Add(new HermitianFermionTerm( new [] { 2, 2 }.ToLadderSequence()), 1.0);

            hamiltonian.Add(new HermitianFermionTerm( new [] { 0, 2 }.ToLadderSequence()), 1.0);
            hamiltonian.Add(new HermitianFermionTerm( new [] { 1, 3 }.ToLadderSequence()), 1.0);
            hamiltonian.Add(new HermitianFermionTerm( new [] { 2, 6 }.ToLadderSequence()), 1.0);

            hamiltonian.Add(new HermitianFermionTerm( new [] { 0, 2, 2, 0 }.ToLadderSequence()), 1.0);
            hamiltonian.Add(new HermitianFermionTerm( new [] { 1, 3, 3, 1 }.ToLadderSequence()), 1.0);
            hamiltonian.Add(new HermitianFermionTerm( new [] { 2, 6, 6, 2 }.ToLadderSequence()), 1.0);

            hamiltonian.Add(new HermitianFermionTerm( new [] { 0, 2, 2, 1 }.ToLadderSequence()), 1.0);
            hamiltonian.Add(new HermitianFermionTerm( new [] { 1, 3, 3, 2 }.ToLadderSequence()), 1.0);
            hamiltonian.Add(new HermitianFermionTerm( new [] { 2, 6, 6, 5 }.ToLadderSequence()), 1.0);

            hamiltonian.Add(new HermitianFermionTerm( new [] { 0, 2, 4, 3 }.ToLadderSequence()), 1.0);
            hamiltonian.Add(new HermitianFermionTerm( new [] { 1, 4, 3, 2 }.ToLadderSequence()), 1.0);
            hamiltonian.Add(new HermitianFermionTerm( new [] { 2, 4, 5, 3 }.ToLadderSequence()), 1.0);
            return hamiltonian;
        }

        [Fact]
        public void JordanWignerOptimizedEvolutionSetTest()
        {
            var fermionHamiltonian = GenerateTestHamiltonian();
            var pauliHamiltonian = fermionHamiltonian.ToPauliHamiltonian(QubitEncoding.JordanWigner);
            var (energyOffset, nSpinOrbitals, termData) = pauliHamiltonian.ToQSharpFormat();
            var (hZ, hZZ, hPQandPQQR, v01234) = termData;
            
            //Assert.Equal(3.0 * 0.5 + 0.25 * 3.0, energyOffset);
            Assert.Contains(new HTerm((new QArray<Int64> { 0 }, new QArray<Double> { -0.5 * 1.0 -0.25 * 1.0 })), termData.Item1, new Convert.HTermArrayComparer());
            Assert.Contains(new HTerm((new QArray<Int64> { 1 }, new QArray<Double> { -0.5 * 1.0 - 0.25 * 1.0 })), termData.Item1, new Convert.HTermArrayComparer());
            Assert.Contains(new HTerm((new QArray<Int64> { 2 }, new QArray<Double> { -0.5 * 1.0 - 2.0 * 0.25 * 1.0 })), termData.Item1, new Convert.HTermArrayComparer());
            Assert.Contains(new HTerm((new QArray<Int64> { 3 }, new QArray<Double> { - 0.25 * 1.0 })), termData.Item1, new Convert.HTermArrayComparer());
            Assert.Contains(new HTerm((new QArray<Int64> { 6 }, new QArray<Double> { - 0.25 * 1.0 })), termData.Item1, new Convert.HTermArrayComparer());

            Assert.Contains(new HTerm((new QArray<Int64> { 0, 2 }, new QArray<Double> { 0.25 * 1.0 })), termData.Item2, new Convert.HTermArrayComparer());
            Assert.Contains(new HTerm((new QArray<Int64> { 1, 3 }, new QArray<Double> { 0.25 * 1.0 })), termData.Item2, new Convert.HTermArrayComparer());
            Assert.Contains(new HTerm((new QArray<Int64> { 2, 6 }, new QArray<Double> { 0.25 * 1.0 })), termData.Item2, new Convert.HTermArrayComparer());

            Assert.Contains(new HTerm((new QArray<Int64> { 0, 2, 2, 1 }, new QArray<Double> { -0.125 * 1.0 })), termData.Item3, new Convert.HTermArrayComparer());
            Assert.Contains(new HTerm((new QArray<Int64> { 1, 3, 3, 2 }, new QArray<Double> { -0.125 * 1.0 })), termData.Item3, new Convert.HTermArrayComparer());
            Assert.Contains(new HTerm((new QArray<Int64> { 2, 6, 6, 5 }, new QArray<Double> { -0.125 * 1.0 })), termData.Item3, new Convert.HTermArrayComparer());

            // Test incomplete
        }
       

    }



}