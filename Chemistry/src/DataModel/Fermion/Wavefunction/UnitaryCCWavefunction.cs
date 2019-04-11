// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using System;
using System.Linq;
using System.Collections.Generic;
using Microsoft.Quantum.Chemistry.LadderOperators;


namespace Microsoft.Quantum.Chemistry.Fermion
{
    // This class is for Fermion terms that are not grouped into Hermitian bunches.
    // Maybe need a stype for quantum state?

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
        public UnitaryCCWavefunction<TNewIndex> ToNewIndex<TNewIndex>(Func<TIndex, TNewIndex> indexFunction)
        where TNewIndex : IEquatable<TNewIndex>, IComparable<TNewIndex>
        => new UnitaryCCWavefunction<TNewIndex>()
        {
            // Be sure to propagate any change in the ladder operators to the coefficient.
            Reference = this.Reference.ToNewIndex(indexFunction),
            Excitations = this.Excitations
                .ToDictionary(kv => new IndexOrderedSequence<TNewIndex>(
                    kv.Key.ToNewIndex(indexFunction).Sequence, 1), kv => kv.Value * (double)kv.Key.Coefficient)
        };
    }

}



