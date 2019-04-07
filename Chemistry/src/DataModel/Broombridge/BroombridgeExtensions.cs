// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using Microsoft.Quantum.Simulation.Core;

using System;
using System.Linq;
using System.Collections.Generic;
using Microsoft.Quantum.Chemistry;
using System.Numerics;

namespace Microsoft.Quantum.Chemistry
{
    /// <summary>
    /// Extensions for converting orbital integrals to fermion terms.
    /// </summary>
    public static partial class Extensions
    {
        /// <summary>
        /// Builds Hamiltonian from Broombridge if data is available.
        /// </summary>
        public static OrbitalIntegralHamiltonian CreateOrbitalIntegralHamiltonian(
            this Broombridge.Current.ProblemDescription broombridge)
        {
            return broombridge.CreateOrbitalIntegralHamiltonian();
        }

        /// <summary>
        /// Builds Hamiltonian from Broombridge orbital integral data.
        /// </summary>
        internal static OrbitalIntegralHamiltonian CreateOrbitalIntegralHamiltonian(
        this Broombridge.DataStructures.HamiltonianData hamiltonianData)
        {
            var hamiltonian = new OrbitalIntegralHamiltonian();
            // This will convert from Broombridge 1-indexing to 0-indexing.
            hamiltonian.AddTerms
                (hamiltonianData.OneElectronIntegrals.Values
                .Select(o => new OrbitalIntegral(o.Item1
                .Select(k => Convert.ToInt32(k - 1)), o.Item2, OrbitalIntegral.Convention.Mulliken)
                .ToCanonicalForm())
                .Distinct());

            // This will convert from Broombridge 1-indexing to 0-indexing.
            // This will convert to Dirac-indexing.
            hamiltonian.AddTerms
                (hamiltonianData.TwoElectronIntegrals.Values
                .Select(o => new OrbitalIntegral(o.Item1
                .Select(k => Convert.ToInt32(k - 1)), o.Item2, OrbitalIntegral.Convention.Mulliken)
                .ToCanonicalForm())
                .Distinct());
            
            return hamiltonian;
        }
    }
}



