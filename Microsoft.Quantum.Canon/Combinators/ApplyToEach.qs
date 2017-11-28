namespace Microsoft.Quantum.Canon {

    /// # See Also
    /// - @"microsoft.quantum.canon.ApplyToEach"
    operation ApplyToEachAC<'T>(singleQubitOperation : ('T => () : Adjoint, Controlled), register : 'T[])  : ()
    {
        body {
            for (idxQubit in 0..(Length(register) - 1)) {
                singleQubitOperation(register[idxQubit]);
            }
        }

        adjoint auto
        controlled auto
        controlled adjoint auto
    }

    /// # See Also
    /// - @"microsoft.quantum.canon.ApplyToEach"
    operation ApplyToEachA<'T>(singleQubitOperation : ('T => ():Adjoint), register : 'T[])  : ()
    {
        body {
            for (idxQubit in 0..(Length(register) - 1)) {
                singleQubitOperation(register[idxQubit]);
            }
        }

        adjoint auto
    }

    /// # See Also
    /// - @"microsoft.quantum.canon.ApplyToEach"
    operation ApplyToEachC<'T>(singleQubitOperation : ('T => ():Controlled), register : 'T[])  : ()
    {
        body {
            for (idxQubit in 0..(Length(register) - 1)) {
                singleQubitOperation(register[idxQubit]);
            }
        }

        controlled auto
    }

    /// # Summary
    /// Applies a single-qubit operation to each qubit in a register.
    ///
    /// # Input
    /// ## singleQubitOperation
    /// Operation to apply to each qubit.
    /// ## register
    /// Array of qubits on which to apply the given operation.
    ///
    /// # Remarks
    /// ## Example
    /// Prepare a three-qubit $\ket{+} state:
    /// ```qsharp
    ///     using (register = Qubit[3]) {
    ///         ApplyToEach(H, register);
    ///     }
    /// ```
    ///
    /// # See Also
    /// - @"microsoft.quantum.canon.ApplyToEachA"
    /// - @"microsoft.quantum.canon.ApplyToEachC"
    /// - @"microsoft.quantum.canon.ApplyToEachAC"
    operation ApplyToEach<'T>(singleQubitOperation : ('T => ()), register : 'T[])  : ()
    {
        body {
            for (idxQubit in 0..(Length(register) - 1)) {
                singleQubitOperation(register[idxQubit]);
            }
        }
    }

}