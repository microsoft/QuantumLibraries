// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using System;
using System.Linq;
using System.Collections.Generic;



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
    }

}



