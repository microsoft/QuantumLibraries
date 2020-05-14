// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using Microsoft.Quantum.Simulation.Core;

using System.Text.RegularExpressions;
using System.Linq;
using System.Collections.Generic;
using System.IO;

using Microsoft.Quantum.Chemistry.OrbitalIntegrals;

namespace Microsoft.Quantum.Chemistry
{
    /// <summary>
    ///      A minimal problem description to be used for interchange
    ///      with legacy formats.
    /// </summary>
    public struct MinimalProblemDescription
    {
        public int NOrbitals { get; set; }
        public int NElectrons { get; set; }
        public double CoulombRepulsion { get; set; }
        public string MiscellaneousInformation { get; set; }
        public OrbitalIntegralHamiltonian OrbitalIntegralHamiltonian { get; set; }

        public Broombridge.V0_1.ArrayQuantity<long, double> IndicesToArrayQuantity(
            TermType.OrbitalIntegral termType,
            string units = "hartree"
        ) => new Broombridge.V0_1.ArrayQuantity<long, double>
        {
            Units = units,
            Format = "sparse",
            Values = OrbitalIntegralHamiltonian
                    .Terms[termType]
                    .Select(
                        termPair => (
                            termPair
                                .Key
                                .OrbitalIndices
                                .Select(idx => (long)(idx + 1))
                                .ToArray(),
                            termPair.Value.Value
                        )
                    )
                    .ToList()
        };
    }
}
