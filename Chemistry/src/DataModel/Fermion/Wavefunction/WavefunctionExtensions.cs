// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using System;
using System.Linq;
using System.Collections.Generic;
using Microsoft.Quantum.Chemistry.OrbitalIntegrals;
using Microsoft.Quantum.Chemistry.Fermion;

namespace Microsoft.Quantum.Chemistry.Fermion
{ 
    public static partial class Extensions
    {

        /// <summary>
        /// Populates the Unitary coupled-cluster wavefunction with all possible spin-preserving exitations
        /// from occupied orbitals to virtual orbitals.
        /// </summary>
        /// <param name="occupiedOrbitals">Occupied orbitals that annihilation operators act on.</param>
        /// <param name="virtualOrbitals">Virtual orbitals that excitation operators act on.</param>
        public void AddAllUCCSDSingletExcitations(List<int> occupiedOrbitals, List<int> virtualOrbitals)
        {
            var singlesSinglet = new[] { (Spin.u, Spin.u), (Spin.d, Spin.d) };
            var doublesSinglet = new[] { (Spin.u, Spin.u, Spin.u, Spin.u), (Spin.u, Spin.d, Spin.d, Spin.u), (Spin.d, Spin.d, Spin.d, Spin.d) };

            // Populate singles excitations.
            foreach (var occupiedIdx in occupiedOrbitals)
            {
                foreach (var virtualIdx in virtualOrbitals)
                {
                    foreach (var spins in singlesSinglet)
                    {
                        Set(new[] { (virtualIdx, spins.Item1), (occupiedIdx, spins.Item2) }.ToLadderSequence(), new Complex());
                    }
                }
            }


            /// <summary>
            /// Set a term of the wavefunction.
            /// </summary>
            /// <param name="term">Index to term to set amplitude of.</param>
            /// <param name="amplitude">Relative amplitude of term.</param>
            public void Set(IndexOrderedSequence<TIndex> term, Complex amplitude)
            {
                Excitations[term] = amplitude * (double)term.Coefficient;
                term.Coefficient = 1;
            }
        }
    }

}



