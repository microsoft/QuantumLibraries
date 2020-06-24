// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using System.Collections.Generic;
using System.Linq;
using Microsoft.Quantum.Chemistry.Generic;

namespace Microsoft.Quantum.Chemistry.Fermion
{
    public partial class FermionHamiltonian : Hamiltonian<TermType.Fermion, HermitianFermionTerm, DoubleCoeff>
    {
        public FermionHamiltonian() : base() { }

        /// <summary>
        /// Method for collecting all distinct system (orbital) indices.
        /// </summary>
        /// <param name="index">Collate orbital indices from this orbital integral.</param>
        public override void AddToSystemIndices(HermitianFermionTerm index)
        {
            foreach (var idx in index.Sequence)
            {
                SystemIndices.Add(idx.Index);
            }
        }
    }
    
}