// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using System;
using System.Linq;

using Microsoft.Quantum.Chemistry.OrbitalIntegrals;

namespace Microsoft.Quantum.Chemistry.Broombridge
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
            this CurrentVersion.ProblemDescription broombridge)
        {
            return broombridge.CreateOrbitalIntegralHamiltonian();
        }

        /// <summary>
        /// Builds Hamiltonian from Broombridge orbital integral data.
        /// </summary>
        internal static OrbitalIntegralHamiltonian CreateOrbitalIntegralHamiltonian(
        this DataStructures.DataStructure.HamiltonianData hamiltonianData)
        {
            var hamiltonian = new OrbitalIntegralHamiltonian();
            // This will convert from Broombridge 1-indexing to 0-indexing.
            hamiltonian.AddTerms
                (hamiltonianData.OneElectronIntegrals.Values
                .Select(o => new OrbitalIntegral(o.Item1
                .Select(k => (int)(k - 1)), o.Item2, OrbitalIntegral.Convention.Mulliken)
                .ToCanonicalForm())
                .Distinct());

            // This will convert from Broombridge 1-indexing to 0-indexing.
            // This will convert to Dirac-indexing.
            hamiltonian.AddTerms
                (hamiltonianData.TwoElectronIntegrals.Values
                .Select(o => new OrbitalIntegral(o.Item1
                .Select(k => (int)(k - 1)), o.Item2, OrbitalIntegral.Convention.Mulliken)
                .ToCanonicalForm())
                .Distinct());
            
            return hamiltonian;
        }
    }
}



