// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {
    open Microsoft.Quantum.Primitive;

    ///<summary>
    /// Applies unitaries that map |0⟩⊗…⊗|0⟩ to |ψ_1⟩⊗…⊗|ψ_n⟩
    /// where |ψ_k⟩ depends on basis[k]. The correspondence between 
    /// value of basis[k] and |ψ_k⟩ is the following:
    /// 0 -- |0⟩
    /// 1 -- |1⟩
    /// 2 -- |+⟩
    /// 3 -- |i⟩ ( +1 eigenstate of Pauli Y )
    ///</summary>
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

    ///<summary>
    /// Checks if the result of applying two operation givenU and expectedU to 
    /// a basis state is the same. The basis state is described by `basis` parameter. 
    /// See FlipToBasis function for more details on this description. 
    ///</summary>
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
            }
        }
    }

    ///<summary> 
    /// Test helper to be used with AssertOperationsEqualInPlace and AssertOperationsEqualReferenced 
    /// if one needs to check that the operation is identity
    ///</summary>
    operation IdentityTestHelper( q : Qubit[] ) : ()
    {
        body {}
        adjoint auto
        controlled auto
        controlled adjoint auto
    }

    ///<summary>
    /// Checks if the operation givenU is equal to the operation givenU on 
    /// the given input size and up to a given tolerance
    ///</summary>
    operation AssertOperationsEqualInPlace( givenU : (Qubit[] => ()), expectedU : (Qubit[] => () : Adjoint ), inputSize : Int ) : ()
    {
        body
        {
            let tolerance = 1e-5;
            let checkOperation = AssertEqualOnBasisVector( _, givenU, expectedU, tolerance );
            IterateThroughCartesianPower(inputSize,4,checkOperation);
        }
    }

    ///<summary>
    /// Checks if the operation givenU is equal to the operation givenU on 
    /// the given input size and up to a given tolerance 
    /// by checking action of the operations only on the vectors from computational basis. 
    /// This is necessary, but not sufficient condition for the equality of two unitaries.
    ///</summary>
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
