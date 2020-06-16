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
            yield return () => Assert.Equal(expected.CoulombRepulsion.Value, actual.CoulombRepulsion.Value);
            yield return () => Assert.Equal(expected.EnergyOffset.Value, actual.EnergyOffset.Value);
            yield return () => Assert.Equal(expected.NElectrons, actual.NElectrons);
            yield return () => Assert.Equal(expected.NOrbitals, actual.NOrbitals);
            yield return () => Assert.Equal(expected.OrbitalIntegralHamiltonian.Terms.Count, actual.OrbitalIntegralHamiltonian.Terms.Count);

            // Check the actual terms themselves.
            // NB: We don't have to worry about actual having expected term kinds, since we
            //     checked Count above.
            foreach (var termType in expected.OrbitalIntegralHamiltonian.Terms.Keys)
            {
                yield return () => Assert.Equal(
                    expected.OrbitalIntegralHamiltonian.Terms[termType].Count,
                    actual.OrbitalIntegralHamiltonian.Terms[termType].Count
                );

                // Now that we've checked that the count is correct, we can check
                // each individual term.
                foreach (var index in expected.OrbitalIntegralHamiltonian.Terms[termType])
                {
                    yield return () => Assert.Equal(
                        index.Value,
                        actual.OrbitalIntegralHamiltonian.Terms[termType][index.Key]
                    );
                }
            }
        }
    }
}
