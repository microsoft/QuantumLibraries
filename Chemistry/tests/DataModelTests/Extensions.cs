// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

#nullable enable

using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

using Microsoft.Quantum.Chemistry.Broombridge;

using Newtonsoft.Json;

using Xunit;
using FluentAssertions;

namespace Microsoft.Quantum.Chemistry.Tests
{
    internal static class Extensions
    {
        internal static void AssertThat(Action what) =>
            what();

        internal static void AssertThat(IEnumerable<Action> whats)
        {
            var exceptions = new List<Exception>();
            foreach (var what in whats)
            {
                try
                {
                    what();
                }
                catch (Exception ex)
                {
                    exceptions.Add(ex);
                }
            }

            if (exceptions.Count == 0)
            {
                return;
            }
            else if (exceptions.Count == 1)
            {
                throw exceptions.Single();
            }
            else
            {
                throw new AggregateException(exceptions);
            }
        }

        internal static Action HasUnits<T>(this Quantity<T> quantity, string units) =>
            () => Assert.Equal(quantity.Units, units);

        internal static IEnumerable<Action> IsEqualTo(this ElectronicStructureProblem expected, ElectronicStructureProblem actual)
        {   
            // TODO: check rest of metadata.
            yield return () =>
                actual.CoulombRepulsion.Value.Should().Be(
                    expected.CoulombRepulsion.Value,
                    because: "Coulomb repulsion should match"
                );
            yield return () => 
                actual.EnergyOffset.Value.Should().Be(
                    expected.EnergyOffset.Value,
                    because: "energy offset should match"
                );
            yield return () =>
                actual.NElectrons.Should().Be(
                    expected.NElectrons,
                    because: "# of electrons should match"
                );
            yield return () =>
                actual.NOrbitals.Should().Be(
                    expected.NOrbitals,
                    because: "# of orbitals should match"
                );
            // yield return () =>
            //     actual.OrbitalIntegralHamiltonian.Terms.Should().ContainKeys(
            //         expected.OrbitalIntegralHamiltonian.Terms.Keys
            //     ).And.HaveCount(
            //         expected.OrbitalIntegralHamiltonian.Terms.Count,
            //         because: "# of kinds of terms should match"
            //     );

            var onlyInActual =
                actual.OrbitalIntegralHamiltonian.Terms.Keys
                .Except(expected.OrbitalIntegralHamiltonian.Terms.Keys);

            // Make sure that any keys that are only in the actual
            // problem are trivial.
            foreach (var termType in onlyInActual)
            {
                var term = actual.OrbitalIntegralHamiltonian.Terms[termType];
                // Make sure that any terms, if present, have coefficient zero.
                yield return () => term.Values.Should().OnlyContain(item => item.Value == 0.0);
            }

            // Check the actual terms themselves.
            foreach (var termType in expected.OrbitalIntegralHamiltonian.Terms.Keys)
            {
                yield return () =>
                    actual.OrbitalIntegralHamiltonian.Terms[termType].Count.Should().Be(
                        expected.OrbitalIntegralHamiltonian.Terms[termType].Count,
                        because: $"# of terms of type {termType} should match"
                    );

                // Now that we've checked that the count is correct, we can check
                // each individual term.
                foreach (var index in expected.OrbitalIntegralHamiltonian.Terms[termType])
                {
                    yield return () =>
                        actual.OrbitalIntegralHamiltonian.Terms[termType][index.Key].Should().Be(
                            index.Value,
                            $"term with index {index.Key} of type {termType} should match"
                        );
                }
            }
        }
    }
}
