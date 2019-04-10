// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using System;
using System.Linq;
using System.Collections.Generic;

using Microsoft.Quantum.Chemistry.Pauli;

namespace Microsoft.Quantum.Chemistry.Fermion
{
    using static Microsoft.Quantum.Chemistry.Extensions;
    /// <summary>
    /// Extensions for converting fermion terms to Pauli terms.
    /// </summary>
    public static partial class Extensions
    {
        /// <summary>
        /// This approximates the Hamiltonian ground state by a greedy algorithm  
        /// that minimizes only the PP term energies. If there are no PP terms,
        /// states will be occupied in lexicographic order.
        /// </summary>
        /// <returns>
        /// Greedy trial state for minimizing Hamiltonian diagonal one-electron energy.
        /// </returns>
        public static SingleCFWavefunction<int> GreedyStatePreparationSCF(this FermionHamiltonian hamiltonian, int nElectrons)
        {
            if (hamiltonian.Terms.ContainsKey(TermType.Fermion.PP))
            {
                var hPPTermSortedByCoeff = hamiltonian.Terms[TermType.Fermion.PP];
                var spinOrbitalIndices = hPPTermSortedByCoeff.OrderBy(o => o.Value).Select(o => o.Key.Sequence.First().Index).Take((int)nElectrons).ToArray();
                return new SingleCFWavefunction<int>(spinOrbitalIndices);
            }
            else
            {
                return new SingleCFWavefunction<int>(Enumerable.Range(0, nElectrons).ToArray());
            }
        }
        
    }
}



