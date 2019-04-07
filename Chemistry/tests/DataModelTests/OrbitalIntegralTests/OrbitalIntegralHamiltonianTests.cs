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
            hamiltonian.AddTerms(orbitalIntegrals.Select(o => (o, o.Coefficient.ToDouble())));
        }

        [Fact]
        void CheckTermPresent()
        {
            var hamiltonian = new OrbitalIntegralHamiltonian();
            var addTerms0 = orbitalIntegrals.Select(o => (o, o.Coefficient.ToDouble())).ToList();
            hamiltonian.AddTerms(addTerms0);

            // Check that all terms present.
            foreach(var term in addTerms0)
            {
                Assert.Equal(term.Item2, hamiltonian.GetTerm(term.o));
            }

            // Now check that indexing is by value and not by reference.
            var newTerms0Copy = new[] {
                new[] { 0,0 },
                new[] { 1,1 },
                new[] { 0,0,0,0 },
                new[] { 0,1,0,1 },
                new[] { 0,1,1,0 },
                new[] { 1,1,1,1 } }.Select(o => new OrbitalIntegral(o));
            foreach (var term in addTerms0.Zip(newTerms0Copy, (a,b) => (a.Item2.Value, b)))
            {
                Assert.Equal(term.Item1, hamiltonian.GetTerm(term.b).Value);
            }


            var orb = new OrbitalIntegral(new[] { 0,1,1,0}, 4.0);
            Assert.Equal(new[] { 0, 1, 1, 0 }, orb.OrbitalIndices);
            
            Assert.Equal(0.663472101.ToDouble(), hamiltonian.terms[orb.GetTermType()][new OrbitalIntegral(orb.OrbitalIndices, orb.Coefficient)]);
        }

        [Fact]
        void CountTerms()
        {
            var hamiltonian = new OrbitalIntegralHamiltonian();
            hamiltonian.AddTerms(orbitalIntegrals.Select(o => (o, o.Coefficient.ToDouble())));
            var oneNorm = hamiltonian.Norm();
            Assert.Equal(orbitalIntegrals.Count(), hamiltonian.CountTerms());

            hamiltonian.AddTerms(orbitalIntegrals.Select(o => (o, o.Coefficient.ToDouble())));
            Assert.Equal(orbitalIntegrals.Count(), hamiltonian.CountTerms());

            hamiltonian.AddTerms(orbitalIntegrals.Select(o => (o, o.Coefficient.ToDouble())));
            Assert.Equal(orbitalIntegrals.Count(), hamiltonian.CountTerms());

            hamiltonian.AddTerms(orbitalIntegrals.Select(o => (o, o.Coefficient.ToDouble())));
            Assert.Equal(orbitalIntegrals.Count(), hamiltonian.CountTerms());
        }

        [Fact]
        void NormTerms()
        {
            var hamiltonian = new OrbitalIntegralHamiltonian();
            hamiltonian.AddTerms(orbitalIntegrals.Select(o => (o, o.Coefficient.ToDouble())));
            var oneNorm = hamiltonian.Norm();

            hamiltonian.AddTerms(orbitalIntegrals.Select(o => (o, o.Coefficient.ToDouble())));
            Assert.Equal(oneNorm * 2.0, hamiltonian.Norm());

            hamiltonian.AddTerms(orbitalIntegrals.Select(o => (o, o.Coefficient.ToDouble())));
            Assert.Equal(oneNorm * 3.0, hamiltonian.Norm());

            hamiltonian.AddTerms(orbitalIntegrals.Select(o => (o, o.Coefficient.ToDouble())));
            Assert.Equal(oneNorm * 4.0, hamiltonian.Norm());
        }

    }
}