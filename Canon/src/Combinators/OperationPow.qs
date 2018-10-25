// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.


namespace Microsoft.Quantum.Canon
{
    
    operation OperationPowImpl<'T> (oracle : ('T => Unit), power : Int, target : 'T) : Unit
    {
        for (idxApplication in 0 .. power - 1)
        {
            oracle(target);
        }
    }
    
    
    operation OperationPowImplC<'T> (oracle : ('T => Unit : Controlled), power : Int, target : 'T) : Unit
    {
        body (...)
        {
            for (idxApplication in 0 .. power - 1)
            {
                oracle(target);
            }
        }
        
        controlled distribute;
    }
    
    
    operation OperationPowImplA<'T> (oracle : ('T => Unit : Adjoint), power : Int, target : 'T) : Unit
    {
        body (...)
        {
            for (idxApplication in 0 .. power - 1)
            {
                oracle(target);
            }
        }
        
        adjoint invert;
    }
    
    
    operation OperationPowImplCA<'T> (oracle : ('T => Unit : Controlled, Adjoint), power : Int, target : 'T) : Unit
    {
        body (...)
        {
            for (idxApplication in 0 .. power - 1)
            {
                oracle(target);
            }
        }
        
        adjoint invert;
        controlled distribute;
        controlled adjoint distribute;
    }
    
    
    /// # Summary
    /// Given an operation representing a gate $U$, returns a new operation
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
    /// Given an operation representing a gate $U$, returns a new operation
    /// $U^m$ for a power $m$.
    /// The modifier 'C' indicates that the operation is controllable.
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
    function OperationPowC<'T> (oracle : ('T => Unit : Controlled), power : Int) : ('T => Unit : Controlled)
    {
        return OperationPowImplC(oracle, power, _);
    }
    
    
    /// # Summary
    /// Given an operation representing a gate $U$, returns a new operation
    /// $U^m$ for a power $m$.
    /// The modifier 'A' indicates that the operation is adjointable.
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
    function OperationPowA<'T> (oracle : ('T => Unit : Adjoint), power : Int) : ('T => Unit : Adjoint)
    {
        return OperationPowImplA(oracle, power, _);
    }
    
    
    /// # Summary
    /// Given an operation representing a gate $U$, returns a new operation
    /// $U^m$ for a power $m$.
    /// The modifier 'CA' indicates that the operation is controllable and adjointable.
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
    function OperationPowCA<'T> (oracle : ('T => Unit : Controlled, Adjoint), power : Int) : ('T => Unit : Controlled, Adjoint)
    {
        return OperationPowImplCA(oracle, power, _);
    }
    
}


