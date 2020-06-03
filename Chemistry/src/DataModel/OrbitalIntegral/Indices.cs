// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

#nullable enable

using System;
using System.Linq;
using System.Collections.Generic;

using Microsoft.Quantum.Chemistry;
using Microsoft.Quantum.Chemistry.LadderOperators;

namespace Microsoft.Quantum.Chemistry.OrbitalIntegrals
{

    internal static class IndexConventionConversions
    {
        internal static int[] ConvertIndices(
            IEnumerable<int> indices,
            OrbitalIntegral.Convention from,
            OrbitalIntegral.Convention to
        ) =>
            (from, to) switch
            {
                (OrbitalIntegral.Convention.Mulliken, OrbitalIntegral.Convention.Dirac) =>
                    ConvertMullikenToDirac(indices),
                (OrbitalIntegral.Convention.Dirac, OrbitalIntegral.Convention.Mulliken) =>
                    ConvertDiracToMulliken(indices),
                _ when from == to => indices.ToArray(),
                _ => throw new ArgumentException($"Conversion from {from} to {to} not currently supported.")
            };

        private static int[] ConvertMullikenToDirac(IEnumerable<int> indices) =>
            indices.Count() switch
            {
                2 => indices.Select(o => o).ToArray(),
                4 => new int[]
                {
                    indices.ElementAt(0),
                    indices.ElementAt(2),
                    indices.ElementAt(3),
                    indices.ElementAt(1)
                },
                _ => throw new System.ArgumentException(
                    $"Got indices [{String.Join(", ", indices)}], but Mulliken convention for not 2 or 4 indices is not defined."
                )
            };

        private static int[] ConvertDiracToMulliken(IEnumerable<int> indices) =>
            indices.Count() switch
            {
                2 => indices.Select(o => o).ToArray(),
                4 => new int[]
                {
                    indices.ElementAt(0),
                    indices.ElementAt(2),
                    indices.ElementAt(3),
                    indices.ElementAt(1)
                },
                _ => throw new System.ArgumentException(
                    $"Got indices [{String.Join(", ", indices)}], but Dirac convention for not 2 or 4 indices is not defined."
                )
            };
    }

}
