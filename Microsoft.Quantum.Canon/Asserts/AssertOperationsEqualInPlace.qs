// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {
    open Microsoft.Quantum.Primitive;

    /// # Summary
    /// Applies unitaries that map $\ket{0}\otimes\cdots\ket{0}$
    /// to $\ket{\psi_0} \otimes \ket{\psi_{n - 1}}$,
    /// where $\ket{\psi_k}$ depends on `basis[k]`.
    ///
    /// The correspondence between 
    /// value of `basis[k]` and $\ket{\psi_k}$ is the following:
    ///
    /// | `basis[k]` | $\ket{\psi_k}$ |
    /// |===|===|
    /// | 0 | $\ket{0}$ |
    /// | 1 | $\ket{1}$ |
    /// | 2 | $\ket{+}$ |
    /// | 3 | $\ket{i}$ ( +1 eigenstate of Pauli Y ) |
    ///
    /// # Input
    /// ## qubits
    /// Qubit to be operated on.
    /// ## basis
    /// Array of single-qubit basis state IDs (0 <= id <= 3), one for each element of
    /// qubits.
    operation FlipToBasis( qubits : Qubit[], basis : Int[] ) : ()
    {
        body
        {
            if( Length(qubits) != Length(basis) )
            {
                fail "qubits and stateIds must have the same length";
            }
            for( i in 0 .. Length(qubits) - 1 )
            {
                let id = basis[i];
                if( id < 0 || id > 3 )
                {
                    fail "Invalid values in the stateIds array. Must be between 0 and 3";
                }

                if( id == 0 )
                {
                    I(qubits[i]);
                }
                elif( id == 1 )
                {
                    X(qubits[i]);
                }
                elif( id == 2 )
                {
                    H(qubits[i]);
                }
                else
                {
                    H(qubits[i]);
                    S(qubits[i]);
                }
            }
        }
        adjoint auto
    }

    /// # Summary
    /// Checks if the result of applying two operations givenU and expectedU to
    /// a basis state is the same. The basis state is described by `basis` parameter. 
    /// See @"Microsoft.Quantum.Canon.FlipToBasis" function for more details on this
    /// description.
    ///
    /// # Input
    /// ## basis
    /// Basis state specified by a single-qubit basis state ID (0 <= id <= 3) for each of
    /// n qubits.
    /// ## givenU
    /// Operation on n qubits to be checked.
    /// ## expectedU
    /// Reference operation on n qubits that givenU is to be compared against.
    /// ## tolerance
    /// Tolerance for the comparison of probabilities between the two operations'
    /// outcomes.
    operation AssertEqualOnBasisVector( basis : Int[] , givenU : (Qubit[] => ()), expectedU : (Qubit[] => () : Adjoint ), tolerance : Double ) : ()
    {
        body
        {
            using( qubits = Qubit[Length(basis)] )
            {
                for( i in 0 .. Length(basis) - 1 )
                {
                    AssertProb([PauliZ],[qubits[i]],Zero,1.0,"Expecting qubits to be all zero at the beginning", tolerance);
                }
                FlipToBasis(qubits,basis);
                givenU(qubits);
                (Adjoint(expectedU))(qubits);
                (Adjoint(FlipToBasis))(qubits, basis);
                for( i in 0 .. Length(basis) - 1 )
                {
                    AssertProb([PauliZ],[qubits[i]],Zero,1.0,"State must be |0⟩ if givenU is equal to expectedU", tolerance);
                }

                ResetAll(qubits);
            }
        }
    }

    /// # Summary
    /// Test helper to be used with AssertOperationsEqualInPlace and AssertOperationsEqualReferenced 
    /// if one needs to check that the operation is identity.
    ///
    /// # Input
    /// ## q
    /// Qubit register to be operated upon.
    operation IdentityTestHelper( q : Qubit[] ) : ()
    {
        body {}
        adjoint auto
        controlled auto
        controlled adjoint auto
    }

    /// # Summary
    /// Checks if the operation givenU is equal to the operation expectedU on
    /// the given input size.
    ///
    /// # Input
    /// ## givenU
    /// Operation on n qubits to be checked.
    /// ## expectedU
    /// Reference operation on n qubits that givenU is to be compared against.
    /// ## inputSize
    /// The number of qubits n that the operations givenU and expectedU operate on.
    operation AssertOperationsEqualInPlace( givenU : (Qubit[] => ()), expectedU : (Qubit[] => () : Adjoint ), inputSize : Int ) : ()
    {
        body
        {
            let tolerance = 1e-5;
            let checkOperation = AssertEqualOnBasisVector( _, givenU, expectedU, tolerance );
            IterateThroughCartesianPower(inputSize,4,checkOperation);
        }
    }

    /// # Summary
    /// Checks if the operation givenU is equal to the operation expectedU on
    /// the given input size  by checking the action of the operations only on
    /// the vectors from the computational basis.
    /// This is a necessary, but not sufficient, condition for the equality of
    /// two unitaries.
    ///
    /// # Input
    /// ## givenU
    /// Operation on n qubits to be checked.
    /// ## expectedU
    /// Reference operation on n qubits that givenU is to be compared against.
    /// ## inputSize
    /// The number of qubits n that the operations givenU and expectedU operate on.
    operation AssertOperationsEqualInPlaceCompBasis( givenU : (Qubit[] => ()), expectedU : (Qubit[] => () : Adjoint ), inputSize : Int ) : ()
    {
        body
        {
            let tolerance = 1e-5;
            let checkOperation = AssertEqualOnBasisVector( _, givenU, expectedU, tolerance );
            IterateThroughCartesianPower(inputSize,2,checkOperation);
        }
    }

}
