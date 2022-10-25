// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using System;
using System.Collections.Generic;

using Microsoft.Quantum.Chemistry.Generic;

namespace Microsoft.Quantum.Chemistry.OrbitalIntegrals
{
 
    public class OrbitalIntegralHamiltonian : Hamiltonian<TermType.OrbitalIntegral, OrbitalIntegral, DoubleCoeff>
    {
        /// <summary>
        /// Constructor for empty orbital integral Hamiltonian.
        /// </summary>
        public OrbitalIntegralHamiltonian() : base() {}
        
        /// <summary>
        /// Constructor for orbital integral Hamiltonian from list of 
        /// orbital integral terms.
        /// </summary>
        /// <param name="terms">Terms to add to Hamiltonian.</param>
        public OrbitalIntegralHamiltonian(IEnumerable<OrbitalIntegral> terms) : base()
        {
            Add(terms);
        }

        /// <summary>
        /// Method for adding an orbital integral term to a Hamiltonian.
        /// </summary>
        /// <param name="orbitalIntegral">Orbital integral to add to Hamiltonian.</param>
        public void Add(OrbitalIntegral orbitalIntegral)
        {
            Add(new OrbitalIntegral(orbitalIntegral.OrbitalIndices, symmetry: orbitalIntegral.Symmetry), orbitalIntegral.Coefficient);
        }

        /// <summary>
        /// Method for adding multiple orbital integrals to a Hamiltonian.
        /// </summary>
        /// <param name="orbitalIntegral">Orbital integrals to add to Hamiltonian.</param>
        public void Add(IEnumerable<OrbitalIntegral> terms)
        {
            foreach (var term in terms)
            {
                Add(term);
            }
        }

        /// <summary>
        /// Method for collecting all distinct system (orbital) indices.
        /// </summary>
        /// <param name="orbitalIntegral">Collate orbital indices from this orbital integral.</param>
        public override void AddToSystemIndices(OrbitalIntegral index)
        {
            foreach(var idx in index.OrbitalIndices)
            {
                SystemIndices.Add(idx);
            }
        }
        
    }

     
}
 