// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using Microsoft.Quantum.Simulation.Core;

using System;
using System.Linq;
using System.Collections.Generic;

using System.Runtime.Serialization.Formatters.Binary;
using System.IO;
using System.IO.Compression;
using YamlDotNet.Serialization;
using Microsoft.Extensions.Logging;

namespace Microsoft.Quantum.Chemistry
{

    public class OrbitalIntegralHamiltonian : Hamiltonian<TermType.OrbitalIntegral, OrbitalIntegral>
    {
        /// <summary>
        /// Constructor for empty orbital integral Hamiltonian.
        /// </summary>
        public OrbitalIntegralHamiltonian() : base() { }

        /// <summary>
        /// Method for adding an orbital integral term to a Hamiltonian.
        /// </summary>
        /// <param name="orbitalIntegral">Orbital integral to add to Hamiltonian.</param>
        public void AddTerm(OrbitalIntegral orbitalIntegral)
        {
            AddTerm(orbitalIntegral, orbitalIntegral.Coefficient);
        }

        /// <summary>
        /// Method for adding multiple orbital integrals to a Hamiltonian.
        /// </summary>
        /// <param name="orbitalIntegral">Orbital integrals to add to Hamiltonian.</param>
        public void AddTerms(IEnumerable<OrbitalIntegral> terms)
        {
            foreach (var term in terms)
            {
                AddTerm(term);
            }
        }

        public override void AddToSystemIndices(OrbitalIntegral index)
        {
            foreach(var idx in index.OrbitalIndices)
            {
                systemIndices.Add(Convert.ToInt32(idx));
            }
        }
    }
}