// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using Xunit;
using Microsoft.Quantum.Chemistry;
using Microsoft.Quantum.Simulation.Core;

using System.Text.RegularExpressions;
using System;
using System.Linq;
using System.Collections.Generic;

namespace Microsoft.Quantum.Chemistry.Tests
{
    using static FermionTermType.Common;
    using FermionTerm = FermionTerm;
    using FermionTermType = FermionTermType;
    using SpinOrbital = SpinOrbital;

    public class JordanWignerOptimizedEvolutionSetTests
    {
        public FermionHamiltonian GenerateTestHamiltonian()
        {

            Int64 nOrbitals = 6;
            Int64 nElectrons = 2;
            Dictionary<FermionTermType, List<FermionTerm>> fermionTerms = new Dictionary<FermionTermType, List<FermionTerm>>
            {
                { PPTermType, new List<FermionTerm>() },
                { PQTermType, new List<FermionTerm>() },
                { PQQPTermType, new List<FermionTerm>() },
                { PQQRTermType, new List<FermionTerm>() },
                { PQRSTermType, new List<FermionTerm>() }
            };

            fermionTerms[PPTermType].Add(new FermionTerm(nOrbitals, new Int64[] { 1, 0 }, new Int64[] { 0, 0 }, 1.0));
            fermionTerms[PPTermType].Add(new FermionTerm(nOrbitals, new Int64[] { 1, 0 }, new Int64[] { 1, 1 }, 1.0));
            fermionTerms[PPTermType].Add(new FermionTerm(nOrbitals, new Int64[] { 1, 0 }, new Int64[] { 2, 2 }, 1.0));

            fermionTerms[PQTermType].Add(new FermionTerm(nOrbitals, new Int64[] { 1, 0 }, new Int64[] { 0, 2 }, 1.0));
            fermionTerms[PQTermType].Add(new FermionTerm(nOrbitals, new Int64[] { 1, 0 }, new Int64[] { 1, 3 }, 1.0));
            fermionTerms[PQTermType].Add(new FermionTerm(nOrbitals, new Int64[] { 1, 0 }, new Int64[] { 2, 6 }, 1.0));

            fermionTerms[PQQPTermType].Add(new FermionTerm(nOrbitals, new Int64[] { 1, 1, 0, 0 }, new Int64[] { 0, 2, 2, 0 }, 1.0));
            fermionTerms[PQQPTermType].Add(new FermionTerm(nOrbitals, new Int64[] { 1, 1, 0, 0 }, new Int64[] { 1, 3, 3, 1 }, 1.0));
            fermionTerms[PQQPTermType].Add(new FermionTerm(nOrbitals, new Int64[] { 1, 1, 0, 0 }, new Int64[] { 2, 6, 6, 2 }, 1.0));

            fermionTerms[PQQRTermType].Add(new FermionTerm(nOrbitals, new Int64[] { 1, 1, 0, 0 }, new Int64[] { 0, 2, 2, 1 }, 1.0));
            fermionTerms[PQQRTermType].Add(new FermionTerm(nOrbitals, new Int64[] { 1, 1, 0, 0 }, new Int64[] { 1, 3, 3, 2 }, 1.0));
            fermionTerms[PQQRTermType].Add(new FermionTerm(nOrbitals, new Int64[] { 1, 1, 0, 0 }, new Int64[] { 2, 6, 6, 5 }, 1.0));

            fermionTerms[PQRSTermType].Add(new FermionTerm(nOrbitals, new Int64[] { 1, 1, 0, 0 }, new Int64[] { 0, 2, 4, 3 }, 1.0));
            fermionTerms[PQRSTermType].Add(new FermionTerm(nOrbitals, new Int64[] { 1, 1, 0, 0 }, new Int64[] { 1, 4, 3, 2 }, 1.0));
            fermionTerms[PQRSTermType].Add(new FermionTerm(nOrbitals, new Int64[] { 1, 1, 0, 0 }, new Int64[] { 2, 4, 5, 3 }, 1.0));
            return new FermionHamiltonian(fermionTerms, nOrbitals, nElectrons: nElectrons);
        }

        [Fact]
        public void JordanWignerOptimizedEvolutionSetTest()
        {
            var hamiltonian = JordanWignerEncoding.Create(GenerateTestHamiltonian());
            var termData = hamiltonian.Terms;
            
            Assert.True(GenerateTestHamiltonian().VerifyFermionTerms());
            Assert.Equal(hamiltonian.energyOffset, 3.0 * 0.5 + 0.25 * 3.0);
            Assert.Contains(new HTerm((new QArray<Int64> { 0 }, new QArray<Double> { -0.5 * 1.0 -0.25 * 1.0 })), termData.Item1, new Extensions.HTermArrayComparer());
            Assert.Contains(new HTerm((new QArray<Int64> { 1 }, new QArray<Double> { -0.5 * 1.0 - 0.25 * 1.0 })), termData.Item1, new Extensions.HTermArrayComparer());
            Assert.Contains(new HTerm((new QArray<Int64> { 2 }, new QArray<Double> { -0.5 * 1.0 - 2.0 * 0.25 * 1.0 })), termData.Item1, new Extensions.HTermArrayComparer());
            Assert.Contains(new HTerm((new QArray<Int64> { 3 }, new QArray<Double> { - 0.25 * 1.0 })), termData.Item1, new Extensions.HTermArrayComparer());
            Assert.Contains(new HTerm((new QArray<Int64> { 6 }, new QArray<Double> { - 0.25 * 1.0 })), termData.Item1, new Extensions.HTermArrayComparer());

            Assert.Contains(new HTerm((new QArray<Int64> { 0, 2 }, new QArray<Double> { 0.25 * 1.0 })), termData.Item2, new Extensions.HTermArrayComparer());
            Assert.Contains(new HTerm((new QArray<Int64> { 1, 3 }, new QArray<Double> { 0.25 * 1.0 })), termData.Item2, new Extensions.HTermArrayComparer());
            Assert.Contains(new HTerm((new QArray<Int64> { 2, 6 }, new QArray<Double> { 0.25 * 1.0 })), termData.Item2, new Extensions.HTermArrayComparer());

            Assert.Contains(new HTerm((new QArray<Int64> { 0, 2, 2, 1 }, new QArray<Double> { -0.125 * 1.0 })), termData.Item3, new Extensions.HTermArrayComparer());
            Assert.Contains(new HTerm((new QArray<Int64> { 1, 3, 3, 2 }, new QArray<Double> { -0.125 * 1.0 })), termData.Item3, new Extensions.HTermArrayComparer());
            Assert.Contains(new HTerm((new QArray<Int64> { 2, 6, 6, 5 }, new QArray<Double> { -0.125 * 1.0 })), termData.Item3, new Extensions.HTermArrayComparer());

            // Test incomplete
        }

        [Fact]
        public void HamiltonianAccumulateTermsTest()
        {
            var hamiltonian0 = new FermionHamiltonian(10, 1);
            hamiltonian0.AddFermionTerm(new OrbitalIntegral(new int[] { 0, 1, 2, 0 }, 1.0));

            var hamiltonian1 = new FermionHamiltonian(10, 1);
            hamiltonian1.AddFermionTerm(new OrbitalIntegral(new int[] { 0, 2, 1, 0 }, 1.0));

            Assert.True(hamiltonian0.FermionTerms[PQQRTermType].SequenceEqual(hamiltonian1.FermionTerms[PQQRTermType]));

            var hamiltonian2 = new FermionHamiltonian(10, 1);
            hamiltonian2.AddFermionTerm(new OrbitalIntegral(new int[] { 3, 1, 2, 3 }, 1.0));

            var hamiltonian3 = new FermionHamiltonian(10, 1);
            hamiltonian3.AddFermionTerm(new OrbitalIntegral(new int[] { 3, 2, 1, 3 }, 1.0));
            Assert.True(hamiltonian2.FermionTerms[PQQRTermType].SequenceEqual(hamiltonian3.FermionTerms[PQQRTermType]));

            var hamiltonian4 = new FermionHamiltonian(10, 1);
            hamiltonian4.AddFermionTerm(new OrbitalIntegral(new int[] { 3, 1, 2, 3 }, 3.0));
            hamiltonian4.SortAndAccumulate();

            var hamiltonian5 = new FermionHamiltonian(10, 1);
            hamiltonian5.AddFermionTerm(new OrbitalIntegral(new int[] { 3, 2, 1, 3 }, 2.0));
            hamiltonian5.AddFermionTerm(new OrbitalIntegral(new int[] { 3, 2, 1, 3 }, 1.0));
            hamiltonian5.SortAndAccumulate();
            Assert.True(hamiltonian4.FermionTerms[PQQRTermType].SequenceEqual(hamiltonian5.FermionTerms[PQQRTermType]));

            var hamiltonian6 = new FermionHamiltonian(10, 1);
            hamiltonian6.AddFermionTerm(new OrbitalIntegral(new int[] { 3, 1, 2, 3 }, 1.0));
            hamiltonian6.AddFermionTerm(new OrbitalIntegral(new int[] { 3, 2, 1, 3 }, 1.0));
            hamiltonian6.AddFermionTerm(new OrbitalIntegral(new int[] { 3, 2, 1, 3 }, 1.0));
            hamiltonian6.SortAndAccumulate();
            Assert.True(hamiltonian4.FermionTerms[PQQRTermType].SequenceEqual(hamiltonian6.FermionTerms[PQQRTermType]));
        }

    }



}