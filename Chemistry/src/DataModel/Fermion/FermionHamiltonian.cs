// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using Microsoft.Quantum.Chemistry.Hamiltonian;

namespace Microsoft.Quantum.Chemistry.Fermion
{
    public class FermionHamiltonian : GenericHamiltonian<TermType.Fermion, FermionTermHermitian, Double>
    {
        public FermionHamiltonian() : base() { }
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