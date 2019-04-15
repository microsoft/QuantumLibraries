// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using System;
using System.Linq;
using System.Collections.Generic;
using Microsoft.Quantum.Chemistry.LadderOperators;
using System.Numerics;
using Newtonsoft.Json;
using Microsoft.Quantum.Chemistry;

namespace Microsoft.Quantum.Chemistry.Fermion
{
    
    /// <summary>
    /// The sparse multi-configurational wavefunction is a superposition of a small number
    /// of single-configurational wavefunctions. In general, the quantum gate complexity of
    /// preparing this state is at least linear in this number.
    /// </summary>
    /// <typeparam name="TIndex">Index of fermion ladder operator.</typeparam>
    public class SparseMultiCFWavefunction<TIndex>
        where TIndex : IEquatable<TIndex>, IComparable<TIndex>
    {
        /// <summary>
        /// Reference state that excitations act on.
        /// </summary>
        [JsonConverter(typeof(Json.LadderSequenceJsonConverter))]
        public SingleCFWavefunction<TIndex> Reference { get; set; }

        /// <summary>
        /// Un-normalized amplitudes of excitations applied to reference state.
        /// </summary>
        [JsonConverter(typeof(Json.FermionWavefunctionJsonConverter))]
        public Dictionary<IndexOrderedSequence<TIndex>, Complex> Excitations { get; set; }


        #region Constructors
        /// <summary>
        /// Constructor for empty instance.
        /// </summary>
        public SparseMultiCFWavefunction() 
        {
            Reference = new SingleCFWavefunction<TIndex>();
            Excitations = new Dictionary<IndexOrderedSequence<TIndex>, Complex>();
        }
        #endregion

        /// <summary>
        /// Changes the indexing scheme of this instance.
        /// </summary>
        /// <typeparam name="TNewIndex">Type of the new indexing scheme.</typeparam>
        /// <param name="indexFunction">Function for mapping the current scheme to the new scheme.</param>
        /// <returns>Instance with a new index type.</returns>
        public SparseMultiCFWavefunction<TNewIndex> ToNewIndex<TNewIndex>(Func<TIndex, TNewIndex> indexFunction)
        where TNewIndex : IEquatable<TNewIndex>, IComparable<TNewIndex>
        => new SparseMultiCFWavefunction<TNewIndex>()
        {
            Reference = this.Reference.ToNewIndex(indexFunction),
            Excitations = this.Excitations
                .ToDictionary(kv => new IndexOrderedSequence<TNewIndex>(
                    kv.Key.ToNewIndex(indexFunction).Sequence, 1), kv => kv.Value * (double)kv.Key.Coefficient)
        };


        /// <summary>
        /// Set a term of the wavefunction.
        /// </summary>
        /// <param name="index">Index to term to set.</param>
        /// <param name="coefficient">Relative amplitude of term.</param>
        public void Set(Complex amplitude, IndexOrderedSequence<TIndex> term)
        {
            Excitations[term] = amplitude * (double) term.Coefficient;
            term.Coefficient = 1;
        }


    }

}



