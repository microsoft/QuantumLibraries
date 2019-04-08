// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using Microsoft.Quantum.Chemistry.Generic;

namespace Microsoft.Quantum.Chemistry.Fermion
{
    public class FermionHamiltonian : Hamiltonian<TermType.Fermion, FermionTermHermitian, DoubleCoeff>
    {
        public FermionHamiltonian() : base() { }

        /// <summary>
        /// Method for collecting all distinct system (orbital) indices.
        /// </summary>
        /// <param name="orbitalIntegral">Collate orbital indices from this orbital integral.</param>
        public override void AddToSystemIndices(FermionTermHermitian index)
        {
            foreach (var idx in index.Sequence)
            {
                SystemIndices.Add(idx.Index);
            }
        }
    }
    
}