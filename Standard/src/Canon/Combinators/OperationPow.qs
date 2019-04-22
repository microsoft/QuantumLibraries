// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon
{
    
    operation OperationPowImpl<'T> (oracle : ('T => Unit), power : Int, target : 'T) : Unit {
        for (idxApplication in 0 .. power - 1)
        {
            oracle(target);
        }
    }
    
    
    operation OperationPowImplC<'T> (oracle : ('T => Unit is Ctl), power : Int, target : 'T) : Unit is Ctl {
        for (idxApplication in 0 .. power - 1)
        {
            oracle(target);
        }
    }
    
    
    operation OperationPowImplA<'T> (oracle : ('T => Unit is Adj), power : Int, target : 'T) : Unit is Adj 
    {
        for (idxApplication in 0 .. power - 1)
        {
            oracle(target);
        }
    }
    
    
    operation OperationPowImplCA<'T> (oracle : ('T => Unit is Adj + Ctl), power : Int, target : 'T) : Unit is Adj + Ctl
    {
        for (idxApplication in 0 .. power - 1)
        {
            oracle(target);
        }
    }
    
    
    /// # Summary
	/// Raises an operation to a power.
	/// 
    /// That is, given an operation representing a gate $U$, returns a new operation
    /// $U^m$ for a power $m$.
    ///
    /// # Input
    /// ## oracle
    /// An operation $U$ representing the gate to be repeated.
    /// ## power
    /// The number of times that $U$ is to be repeated.
    ///
    /// # Output
    /// A new operation representing $U^m$, where $m = \texttt{power}$.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The type of the operation to be powered.
    ///
    /// # See Also
    /// - @"microsoft.quantum.canon.operationpowc"
    /// - @"microsoft.quantum.canon.operationpowa"
    /// - @"microsoft.quantum.canon.operationpowca"
    function OperationPow<'T> (oracle : ('T => Unit), power : Int) : ('T => Unit)
    {
        return OperationPowImpl(oracle, power, _);
    }
    
    
    /// # Summary
	/// Raises an operation to a power.
    /// The modifier `C` indicates that the operation is controllable.
	/// 
    /// That is, given an operation representing a gate $U$, returns a new operation
    /// $U^m$ for a power $m$.
    ///
    /// # Input
    /// ## oracle
    /// An operation $U$ representing the gate to be repeated.
    /// ## power
    /// The number of times that $U$ is to be repeated.
    ///
    /// # Output
    /// A new operation representing $U^m$, where $m = \texttt{power}$.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The type of the operation to be powered.
    ///
    /// # See Also
    /// - @"microsoft.quantum.canon.operationpow"
    function OperationPowC<'T> (oracle : ('T => Unit is Ctl), power : Int) : ('T => Unit is Ctl)
    {
        return OperationPowImplC(oracle, power, _);
    }
    
    
    /// # Summary
	/// Raises an operation to a power.
    /// The modifier `A` indicates that the operation is adjointable.
	/// 
    /// That is, given an operation representing a gate $U$, returns a new operation
    /// $U^m$ for a power $m$.
    ///
    /// # Input
    /// ## oracle
    /// An operation $U$ representing the gate to be repeated.
    /// ## power
    /// The number of times that $U$ is to be repeated.
    ///
    /// # Output
    /// A new operation representing $U^m$, where $m = \texttt{power}$.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The type of the operation to be powered.
    ///
    /// # See Also
    /// - @"microsoft.quantum.canon.operationpow"
    function OperationPowA<'T> (oracle : ('T => Unit is Adj), power : Int) : ('T => Unit is Adj)
    {
        return OperationPowImplA(oracle, power, _);
    }
    
    
    /// # Summary
	/// Raises an operation to a power.
    /// The modifier `A` indicates that the operation is controllable and adjointable.
	/// 
    /// That is, given an operation representing a gate $U$, returns a new operation
    /// $U^m$ for a power $m$.
    ///
    /// # Input
    /// ## oracle
    /// An operation $U$ representing the gate to be repeated.
    /// ## power
    /// The number of times that $U$ is to be repeated.
    ///
    /// # Output
    /// A new operation representing $U^m$, where $m = \texttt{power}$.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The type of the operation to be powered.
    ///
    /// # See Also
    /// - @"microsoft.quantum.canon.operationpow"
    function OperationPowCA<'T> (oracle : ('T => Unit is Ctl + Adj), power : Int) : ('T => Unit is Ctl + Adj)
    {
        return OperationPowImplCA(oracle, power, _);
    }
    
}


