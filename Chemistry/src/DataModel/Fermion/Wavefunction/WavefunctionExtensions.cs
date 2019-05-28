// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using System;
using System.Linq;
using System.Numerics;
using System.Collections.Generic;
using Microsoft.Quantum.Chemistry.OrbitalIntegrals;
using Microsoft.Quantum.Chemistry.Fermion;
using Microsoft.Quantum.Chemistry.LadderOperators;

namespace Microsoft.Quantum.Chemistry.Fermion
{ 
    public static partial class Extensions
    {
        /// <summary>
        /// Populates the Unitary coupled-cluster wavefunction with all possible spin-preserving exitations
        /// from occupied orbitals to virtual orbitals.
        /// </summary>
        /// <param name="occupiedOrbitals">Occupied orbitals that annihilation operators act on.</param>
        /// <param name="nOrbitals">Number of orbitals.</param>
        public static UnitaryCCWavefunction<SpinOrbital> AddAllUCCSDSingletExcitations(this SingleCFWavefunction<SpinOrbital> occupiedOrbitals, int nOrbitals)
        {
            // Create list of unoccupied spin-orbitals
            var virtualOrbitals = Enumerable.Range(0, nOrbitals).Select(x => new SpinOrbital(x, Spin.u))
                .Concat(Enumerable.Range(0, nOrbitals).Select(x => new SpinOrbital(x, Spin.d)))
                .Except(occupiedOrbitals.Sequence.Select(o => o.Index));

            var outputState = new UnitaryCCWavefunction<SpinOrbital>()
            {
                Reference = occupiedOrbitals
            };

            // Populate singles excitations.
            foreach (var occupiedIdx in occupiedOrbitals.Sequence.Select(o => o.Index))
            {
                foreach (var virtualIdx in virtualOrbitals)
                {
                    if (occupiedIdx.Spin == virtualIdx.Spin)
                    {
                        var term = new IndexOrderedSequence<SpinOrbital>(new[] { virtualIdx, occupiedIdx }.ToLadderSequence());
                        outputState.Set(term, new Complex());
                    }
                }
            }

            // Populate doubles excitations.
            foreach (var occupiedIdx0 in occupiedOrbitals.Sequence.Select(o => o.Index))
            {
                foreach (var occupiedIdx1 in occupiedOrbitals.Sequence.Select(o => o.Index))
                {
                    // Only loop over distinct pairs of occupied spin orbitals
                    if (occupiedIdx0.ToInt() < occupiedIdx1.ToInt())
                    {
                        foreach (var virtualIdx0 in virtualOrbitals)
                        {
                            foreach (var virtualIdx1 in virtualOrbitals)
                            {
                                // Only loop over distinct pairs of occupied spin orbitals
                                if (virtualIdx0.ToInt() < virtualIdx1.ToInt())
                                {
                                    if (virtualIdx0.Spin + virtualIdx1.Spin == occupiedIdx0.Spin + occupiedIdx1.Spin)
                                    {
                                        var term = new IndexOrderedSequence<SpinOrbital>(new[] { virtualIdx1, virtualIdx0, occupiedIdx1, occupiedIdx0 }.ToLadderSequence());
                                        outputState.Set(term, new Complex());
                                    }
                                }
                            }
                        }
                    }
                }
            }
            return outputState;
        }
    }

}



