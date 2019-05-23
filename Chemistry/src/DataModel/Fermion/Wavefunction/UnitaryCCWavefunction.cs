// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using System;
using System.Numerics;
using System.Linq;
using System.Collections.Generic;
using Microsoft.Quantum.Chemistry.LadderOperators;
using Microsoft.Quantum.Chemistry.OrbitalIntegrals;

namespace Microsoft.Quantum.Chemistry.Fermion
{
    // For now, UCC is a subclass of MCF. It should eventually be a Hamiltonian
    // + a WavefunctionSCF.
    // 
    public class UnitaryCCWavefunction<TIndex> : SparseMultiCFWavefunction<TIndex>
        where TIndex : IEquatable<TIndex>, IComparable<TIndex>
    {
        public UnitaryCCWavefunction() : base() { }

        /// <summary>
        /// Changes the indexing scheme of this instance.
        /// </summary>
        /// <typeparam name="TNewIndex">Type of the new indexing scheme.</typeparam>
        /// <param name="indexFunction">Function for mapping the current scheme to the new scheme.</param>
        /// <returns>Instance with a new index type.</returns>
        public UnitaryCCWavefunction<TNewIndex> SelectIndex<TNewIndex>(Func<TIndex, TNewIndex> indexFunction)
        where TNewIndex : IEquatable<TNewIndex>, IComparable<TNewIndex>
        => new UnitaryCCWavefunction<TNewIndex>()
        {
            // Be sure to propagate any change in the ladder operators to the coefficient.
            Reference = this.Reference.SelectIndex(indexFunction),
            Excitations = this.Excitations
                .ToDictionary(kv => new IndexOrderedSequence<TNewIndex>(
                    kv.Key.SelectIndex(indexFunction).Sequence, 1), kv => kv.Value * (double)kv.Key.Coefficient)
        };

        // Create excitations from data structure

        /// <summary>
        /// Populates the Unitary coupled-cluster wavefunction with all possible spin-preserving exitations
        /// from occupied orbitals to virtual orbitals.
        /// </summary>
        /// <param name="occupiedOrbitals">Occupied orbitals that annihilation operators act on.</param>
        /// <param name="virtualOrbitals">Virtual orbitals that excitation operators act on.</param>
        public static UnitaryCCWavefunction<SpinOrbital> AddAllUCCSDSingletExcitations(SingleCFWavefunction<SpinOrbital> occupiedOrbitals, int nOrbitals)
        {

            var singlesSinglet = new[] { (Spin.u, Spin.u), (Spin.d, Spin.d) };
            var doublesSinglet = new[] { (Spin.u, Spin.u, Spin.u, Spin.u), (Spin.u, Spin.d, Spin.d, Spin.u), (Spin.d, Spin.d, Spin.d, Spin.d) };

            var virtualOrbitals = Enumerable.Range(0, nOrbitals).Select(x => new SpinOrbital(x, Spin.u))
                .Concat(Enumerable.Range(0, nOrbitals).Select(x => new SpinOrbital(x, Spin.d)))
                .Except(occupiedOrbitals.Sequence.Select(o => o.Index));
            ;

            var outputState = new UnitaryCCWavefunction<SpinOrbital>()
            {
                Reference = occupiedOrbitals
            };

            // Populate singles excitations.
            foreach (var occupiedIdx in occupiedOrbitals.Sequence.Select(o => o.Index))
            {
                foreach(var virtualIdx in virtualOrbitals)
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



