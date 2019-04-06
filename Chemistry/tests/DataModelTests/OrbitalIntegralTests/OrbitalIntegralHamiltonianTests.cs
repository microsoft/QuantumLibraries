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

    public class OrbitalIntegralHamiltonianTests
    {
        // These orbital integrals represent terms in molecular Hydrogen
        private static IEnumerable<OrbitalIntegral> orbitalIntegrals = new OrbitalIntegral[]
            {
                    new OrbitalIntegral(new[] { 0,0 }, -1.252477495),
                    new OrbitalIntegral(new[] { 1,1 }, -0.475934275),
                    new OrbitalIntegral(new[] { 0,0,0,0 }, 0.674493166),
                    new OrbitalIntegral(new[] { 0,1,0,1 }, 0.181287518),
                    new OrbitalIntegral(new[] { 0,1,1,0 }, 0.663472101),
                    new OrbitalIntegral(new[] { 1,1,1,1 }, 0.697398010)
            };

        [Fact]
        void BuildHamiltonian()
        {
            var hamiltonian = new OrbitalIntegralHamiltonian();
            hamiltonian.AddTerms(orbitalIntegrals.Select(o => (o, o.Coefficient)));
        }

        [Fact]
        void CountTerms()
        {
            var hamiltonian = new OrbitalIntegralHamiltonian();
            hamiltonian.AddTerms(orbitalIntegrals.Select(o => (o, o.Coefficient)));
            var oneNorm = hamiltonian.Norm();
            Assert.Equal(orbitalIntegrals.Count(), hamiltonian.CountTerms());

            hamiltonian.AddTerms(orbitalIntegrals.Select(o => (o, o.Coefficient)));
            Assert.Equal(orbitalIntegrals.Count(), hamiltonian.CountTerms());

            hamiltonian.AddTerms(orbitalIntegrals.Select(o => (o, o.Coefficient)));
            Assert.Equal(orbitalIntegrals.Count(), hamiltonian.CountTerms());

            hamiltonian.AddTerms(orbitalIntegrals.Select(o => (o, o.Coefficient)));
            Assert.Equal(orbitalIntegrals.Count(), hamiltonian.CountTerms());
        }

        [Fact]
        void NormTerms()
        {
            var hamiltonian = new OrbitalIntegralHamiltonian();
            hamiltonian.AddTerms(orbitalIntegrals.Select(o => (o, o.Coefficient)));
            var oneNorm = hamiltonian.Norm();

            hamiltonian.AddTerms(orbitalIntegrals.Select(o => (o, o.Coefficient)));
            Assert.Equal(oneNorm * 2.0, hamiltonian.Norm(), 5);

            hamiltonian.AddTerms(orbitalIntegrals.Select(o => (o, o.Coefficient)));
            Assert.Equal(oneNorm * 3.0, hamiltonian.Norm(), 5);

            hamiltonian.AddTerms(orbitalIntegrals.Select(o => (o, o.Coefficient)));
            Assert.Equal(oneNorm * 4.0, hamiltonian.Norm(), 5);
        }

    }
}