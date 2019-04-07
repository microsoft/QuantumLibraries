// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using Microsoft.Quantum.Chemistry.Hamiltonian;

namespace Microsoft.Quantum.Chemistry.Fermion
{
    public class FermionHamiltonian : GenericHamiltonian<TermType.Fermion, FermionTermHermitian, DoubleCoeff>
    {
        public FermionHamiltonian() : base() { }

        /// <summary>
        /// Method for adding a term to a fermion Hamiltonian.
        /// </summary>
        /// <param name="term">Term to be added.</param>
        /// <param name="coefficient">Coefficient of term.</param>
        public void AddTerm(FermionTermHermitian term, double coefficient)
        {
            AddTerm(term, coefficient.ToDouble());
        }

        /// <summary>
        /// Method for collecting all distinct system (orbital) indices.
        /// </summary>
        /// <param name="orbitalIntegral">Collate orbital indices from this orbital integral.</param>
        public override void AddToSystemIndices(FermionTermHermitian index)
        {
            foreach (var idx in index.sequence)
            {
                systemIndices.Add(idx.index);
            }
        }
    }
    /*
    /// <summary>
    /// Method for adding fermion term to a fermion Hamiltonian using double instead of Double type for the coefficients.
    /// </summary>
    /// <param name="orbitalIntegral">Orbital integral to add to Hamiltonian.</param>
    public void AddTerm(FermionTermHermitian term, double value)
    {
        AddTerm(term, value.ToDouble());
    }

    /// <summary>
    /// Method for adding fermion terms to a Hamiltonian using double instead of Double type for the coefficients.
    /// </summary>
    /// <param name="orbitalIntegral">Orbital integrals to add to Hamiltonian.</param>
    public void AddTerms(IEnumerable<OrbitalIntegral> terms)
    {
        foreach (var term in terms)
        {
            AddTerm(term);
        }
    }*/

}