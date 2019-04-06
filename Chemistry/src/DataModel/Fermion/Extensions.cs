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
    /// Class containing methods for converting orbital integrals to fermion terms.
    /// </summary>
    public static partial class Extensions
    {
        /*
        #region Function for state preparation
        /// <summary>
        /// This approximates the Hamiltonian ground state by a greedy algorithm  
        /// that minimizes only the PP term energies. If there are no PP terms,
        /// states will be occupied in lexicographic order.
        /// </summary>
        /// <returns>
        /// Greedy trial state for minimizing Hamiltonian diagonal one-electron energy.
        /// </returns>
        public InputState GreedyStatePreparation(FermionHamiltonian hamiltonian, int nElectrons)
        {
            var label = "Greedy";
            (Double, Double) coeff = (1.0, 0.0);
            var conjugate = Enumerable.Range(0, nElectrons).Select(o => 1).ToArray();

            if (hamiltonian.terms.ContainsKey(TermType.Fermion.PP))
            {
                var hPPTermSortedByCoeff = FermionTerms[PPTermType];
                var spinOrbitals = hPPTermSortedByCoeff.OrderBy(o => o.coeff).Select(o => o.SpinOrbitalIndices.First()).Take((int)NElectrons).ToArray();
                var fermionTerm = new FermionTerm(conjugate, spinOrbitals, 1.0);
                fermionTerm.ToSpinOrbitalCanonicalOrder();
                fermionTerm.coeff = 1.0;

                var superposition = new ((Double, Double), FermionTerm)[] {
                            (coeff, fermionTerm)
                            };
                return new InputState { type = StateType.Sparse_Multi_Configurational, Label = label, Superposition = superposition };
            }
            else
            {
                var fermionTerm = new FermionTerm(
                    (int)NElectrons,
                    conjugate,
                    Enumerable.Range(0, (int)NElectrons).Select(o => (Int64)o).ToArray(),
                    1.0, IndexConvention);
                var superposition = new ((Double, Double), FermionTerm)[] {
                            (coeff, fermionTerm)
                            };
                return new InputState { type = StateType.Sparse_Multi_Configurational, Label = label, Superposition = superposition };
            }
        }
        #endregion
        */

        #region
        #endregion

    }
    
}



