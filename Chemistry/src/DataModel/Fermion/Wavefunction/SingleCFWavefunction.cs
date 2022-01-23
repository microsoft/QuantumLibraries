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

    /// <summary>
    /// Class representing a sequence of fermionic raising and lowering operators, subject to the additional constraints: 
    /// 1) Normal-ordered, where all raising operators are to the left of all lowering operators.
    /// 2) Index-ordered, where are raising(lowering) operators are in ascending(descending) order.
    /// 3) Contains only creation operators.
    /// </summary>
    public class SingleCFWavefunction<TIndex> : IndexOrderedSequence<TIndex>
        where TIndex : IEquatable<TIndex>, IComparable<TIndex>
    {
        #region Constructors
        /// <summary>
        /// Constructor for empty instance.
        /// </summary>
        public SingleCFWavefunction() : base() { }

        /// <summary>
        /// Construct a copy of the input instance.
        /// </summary>
        /// <param name="term">Sequence of ladder operators.</param>
        internal SingleCFWavefunction(SingleCFWavefunction<TIndex> term) : base(term) { }

        /*
        /// <summary>
        /// Construct instance from a normal-ordered sequence of ladder operators.
        /// </summary>
        /// <param name="ladderOperators">Normal-ordered sequence of ladder operators.</param>
        public WavefunctionFermionSCF(LadderSequence ladderOperators) : base(ladderOperators) {
            ThrowExceptionIfNotOnlyRaising();
        }
        */

        /// <summary>
        /// Construct <see cref="WavefunctionFermionSCF"/> from a sequence of integers.
        /// </summary>
        /// <param name="setSequence">Sequence of integers.</param>
        /// <returns>
        /// Sequence of ladder operators with only creation terms in ascending order.
        /// that are normal-ordered.
        /// </returns>
        public SingleCFWavefunction(IEnumerable<TIndex> indices) : base(indices
            .Select(o => new LadderOperator<TIndex>(RaisingLowering.u, o)))
        {
            ThrowExceptionIfNotOnlyRaising();
        }

        #endregion

        /// <summary>
        /// Changes the indexing scheme of this instance.
        /// </summary>
        /// <typeparam name="TNewIndex">Type of the new indexing scheme.</typeparam>
        /// <param name="indexFunction">Function for mapping the current scheme to the new scheme.</param>
        /// <returns>Instance with a new index type.</returns>
        public SingleCFWavefunction<TNewIndex> SelectIndex<TNewIndex>(Func<TIndex, TNewIndex> indexFunction)
        where TNewIndex : IEquatable<TNewIndex>, IComparable<TNewIndex>
            => new SingleCFWavefunction<TNewIndex>(base.SelectIndex(indexFunction).ToIndices());

        /// <summary>
        /// This throws an ArgumentException if the operators in NormalOrderedLadderSequence are not normal-ordered.
        /// </summary>
        private void ThrowExceptionIfNotOnlyRaising()
        {
            if (Sequence.Where(o => o.Type == RaisingLowering.d).Count() > 0)
            {
                throw new ArgumentException("WavefunctionFermionSCF must contain only raising operators.");
            }
        }
    }
}



