// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon
{
    
    /// # Summary
    /// Applies a pair of operations to a given partition of a register into two parts.
    ///
    /// # Input
    /// ## op
    /// The pair of operations to be applied to the given partition.
    /// ## numberOfQubitsToFirstArgument
    /// Number of qubits from target to put into the first part of the partition.
    /// The remaining qubits constitute the second part of the partition.
    /// ## target
    /// A register of qubits that are being partitioned and operated on by the
    /// given two operation.
    ///
    /// # See Also
    /// - @"microsoft.quantum.canon.applytopartitiona"
    /// - @"microsoft.quantum.canon.applytopartitionc"
    /// - @"microsoft.quantum.canon.applytopartitionca"
    operation ApplyToPartition (op : ((Qubit[], Qubit[]) => Unit), numberOfQubitsToFirstArgument : Int, target : Qubit[]) : Unit
    {
        AssertBoolEqual(numberOfQubitsToFirstArgument >= 0, true, $"numberOfQubitsToFirstArgument must be non-negative");
        AssertBoolEqual(Length(target) >= numberOfQubitsToFirstArgument, true, $"Length(target) must greater or equal to numberOfQubitsToFirstArgument");
        op(target[0 .. numberOfQubitsToFirstArgument - 1], target[numberOfQubitsToFirstArgument .. Length(target) - 1]);
    }
    
    
    /// # Summary
    /// Applies a pair of operations to a given partition of a register into two parts.
    /// The modifier 'A' indicates that the operation is adjointable.
    ///
    /// # Input
    /// ## op
    /// The pair of operations to be applied to the given partition.
    /// ## numberOfQubitsToFirstArgument
    /// Number of qubits from target to put into the first part of the partition.
    /// The remaining qubits constitute the second part of the partition.
    /// ## target
    /// A register of qubits that are being partitioned and operated on by the
    /// given two operation.
    ///
    /// # See Also
    /// - @"microsoft.quantum.canon.applytopartition"
    operation ApplyToPartitionA (op : ((Qubit[], Qubit[]) => Unit : Adjoint), numberOfQubitsToFirstArgument : Int, target : Qubit[]) : Unit
    {
        body (...)
        {
            AssertBoolEqual(numberOfQubitsToFirstArgument >= 0, true, $"numberOfQubitsToFirstArgument must be non-negative");
            AssertBoolEqual(Length(target) >= numberOfQubitsToFirstArgument, true, $"Length(target) must greater or equal to numberOfQubitsToFirstArgument");
            op(target[0 .. numberOfQubitsToFirstArgument - 1], target[numberOfQubitsToFirstArgument .. Length(target) - 1]);
        }
        
        adjoint invert;
    }
    
    
    /// # Summary
    /// Applies a pair of operations to a given partition of a register into two parts.
    /// The modifier 'C' indicates that the operation is controllable.
    ///
    /// # Input
    /// ## op
    /// The pair of operations to be applied to the given partition.
    /// ## numberOfQubitsToFirstArgument
    /// Number of qubits from target to put into the first part of the partition.
    /// The remaining qubits constitute the second part of the partition.
    /// ## target
    /// A register of qubits that are being partitioned and operated on by the
    /// given two operation.
    ///
    /// # See Also
    /// - @"microsoft.quantum.canon.applytopartition"
    operation ApplyToPartitionC (op : ((Qubit[], Qubit[]) => Unit : Controlled), numberOfQubitsToFirstArgument : Int, target : Qubit[]) : Unit
    {
        body (...)
        {
            AssertBoolEqual(numberOfQubitsToFirstArgument >= 0, true, $"numberOfQubitsToFirstArgument must be non-negative");
            AssertBoolEqual(Length(target) >= numberOfQubitsToFirstArgument, true, $"Length(target) must greater or equal to numberOfQubitsToFirstArgument");
            op(target[0 .. numberOfQubitsToFirstArgument - 1], target[numberOfQubitsToFirstArgument .. Length(target) - 1]);
        }
        
        controlled distribute;
    }
    
    
    /// # Summary
    /// Applies a pair of operations to a given partition of a register into two parts.
    /// The modifier 'CA' indicates that the operation is controllable and adjointable.
    ///
    /// # Input
    /// ## op
    /// The pair of operations to be applied to the given partition.
    /// ## numberOfQubitsToFirstArgument
    /// Number of qubits from target to put into the first part of the partition.
    /// The remaining qubits constitute the second part of the partition.
    /// ## target
    /// A register of qubits that are being partitioned and operated on by the
    /// given two operation.
    ///
    /// # See Also
    /// - @"microsoft.quantum.canon.applytopartition"
    operation ApplyToPartitionCA (op : ((Qubit[], Qubit[]) => Unit : Controlled, Adjoint), numberOfQubitsToFirstArgument : Int, target : Qubit[]) : Unit
    {
        body (...)
        {
            AssertBoolEqual(numberOfQubitsToFirstArgument >= 0, true, $"numberOfQubitsToFirstArgument must be non-negative");
            AssertBoolEqual(Length(target) >= numberOfQubitsToFirstArgument, true, $"Length(target) must greater or equal to numberOfQubitsToFirstArgument");
            op(target[0 .. numberOfQubitsToFirstArgument - 1], target[numberOfQubitsToFirstArgument .. Length(target) - 1]);
        }
        
        adjoint invert;
        controlled distribute;
        controlled adjoint distribute;
    }
    
}


