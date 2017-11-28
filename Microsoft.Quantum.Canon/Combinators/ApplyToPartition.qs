// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon
{
    /// # Summary 
    /// Applies an operation to a given partition of the register into two parts
    /// 
    /// # Inputs 
    /// ## op
    /// The operation to apply to the partition
    /// ## numberOfQubitsToFirstArgument
    /// Number of qubits from target to put into the first part of the partition.
    /// The rest go into the second part.
    /// ## target
    /// Qubits that are being partitioned 
    operation ApplyToPartition( 
        op : ( (Qubit[],Qubit[]) => () ),
        numberOfQubitsToFirstArgument : Int,
        target : Qubit[]
        ) : () {
        body{ 
            AssertBoolEqual( numberOfQubitsToFirstArgument >= 0 , true,
                "numberOfQubitsToFirstArgument must be non-negative" );
            AssertBoolEqual( Length(target) >=  numberOfQubitsToFirstArgument, true,
                "Length(target) must greater or equal to numberOfQubitsToFirstArgument" );

            op(
                target[ 0 .. numberOfQubitsToFirstArgument - 1 ],
                target[ numberOfQubitsToFirstArgument .. Length(target) - 1 ] );
        }
    }

    /// # See Also 
    /// @"Microsoft.Quantum.Canon.ApplyToPartition"
    operation ApplyToPartitionA( 
        op : ( (Qubit[],Qubit[]) => () : Adjoint ),
        numberOfQubitsToFirstArgument : Int,
        target : Qubit[]
        ) : () {
        body{ 
            AssertBoolEqual( numberOfQubitsToFirstArgument >= 0 , true,
                "numberOfQubitsToFirstArgument must be non-negative" );
            AssertBoolEqual( Length(target) >=  numberOfQubitsToFirstArgument, true,
                "Length(target) must greater or equal to numberOfQubitsToFirstArgument" );

            op(
                target[ 0 .. numberOfQubitsToFirstArgument - 1 ],
                target[ numberOfQubitsToFirstArgument .. Length(target) - 1 ] );
        }
        adjoint auto
    }
    
    /// # See Also 
    /// @"Microsoft.Quantum.Canon.ApplyToPartition"
    operation ApplyToPartitionC( 
        op : ( (Qubit[],Qubit[]) => () : Controlled ),
        numberOfQubitsToFirstArgument : Int,
        target : Qubit[]
        ) : () {
        body{ 
            AssertBoolEqual( numberOfQubitsToFirstArgument >= 0 , true,
                "numberOfQubitsToFirstArgument must be non-negative" );
            AssertBoolEqual( Length(target) >=  numberOfQubitsToFirstArgument, true,
                "Length(target) must greater or equal to numberOfQubitsToFirstArgument" );

            op(
                target[ 0 .. numberOfQubitsToFirstArgument - 1 ],
                target[ numberOfQubitsToFirstArgument .. Length(target) - 1 ] );
        }
        controlled auto
    }

    /// # See Also 
    /// @"Microsoft.Quantum.Canon.ApplyToPartition"
    operation ApplyToPartitionCA( 
        op : ( (Qubit[],Qubit[]) => () : Controlled, Adjoint ),
        numberOfQubitsToFirstArgument : Int,
        target : Qubit[]
        ) : () {
        body{ 
            AssertBoolEqual( numberOfQubitsToFirstArgument >= 0 , true,
                "numberOfQubitsToFirstArgument must be non-negative" );
            AssertBoolEqual( Length(target) >=  numberOfQubitsToFirstArgument, true,
                "Length(target) must greater or equal to numberOfQubitsToFirstArgument" );

            op(
                target[ 0 .. numberOfQubitsToFirstArgument - 1 ],
                target[ numberOfQubitsToFirstArgument .. Length(target) - 1 ] );
        }
        adjoint auto
        controlled auto
        controlled adjoint auto
    }
}
