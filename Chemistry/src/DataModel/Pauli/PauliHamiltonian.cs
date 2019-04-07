// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using Microsoft.Quantum.Chemistry.Hamiltonian;

namespace Microsoft.Quantum.Chemistry.Pauli
{

    /// <summary>
    /// <para>
    /// Jordan-Wigner representation of a general Fermion Hamiltonian <see cref="FermionHamiltonian"/>.
    /// This representation may only be created from instances of <see cref="FermionHamiltonian"/>,
    /// and stores term data in a format suitable for consumption by Q#,
    /// and optimized for a product formula Hamiltonian simulation algorithm.
    /// </para>
    /// <para>
    /// This supports the following <see cref="FermionTermType"/>: 
    /// <see cref="IdentityTermType"/>,
    /// <see cref="PPTermType"/>,
    /// <see cref="PQTermType"/>,
    /// <see cref="PQQPTermType"/>,
    /// <see cref="PQQRTermType"/>,
    /// <see cref="PQRSTermType"/>.
    /// </para>
    /// <para>
    /// Some optimizations are performed: 
    /// <list type="bullet">
    /// <item>PQQP and PP terms are merged where needed.</item>
    /// <item>PQQR and PR terms are merged where possible.</item>
    /// <item>
    /// All PQRS terms with the same set of spin-orbital indices are performed simultaneously,
    /// and only the XXXX, XXYY, XYXY, YXXY, YYYY, YYXX, YXYX, XYYX terms are only performed as
    /// needed.
    /// </item>
    /// <item>Terms in each group of term types are applied in lexicographic ordering.</item>
    /// </list>
    /// </para>
    /// </summary>
    public class PauliHamiltonian : GenericHamiltonian<TermType.PauliTerm, PauliTerm, PauliTermValue>
    {
        public PauliHamiltonian() : base() { }
    }
}
 