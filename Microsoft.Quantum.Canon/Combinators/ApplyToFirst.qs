namespace Microsoft.Quantum.Canon {
    
    // # Overview 
    // Variants of the operation that applies given one, two and three qubit 
    // operation to the first one, two and three qubits of a register

    /// # Summary
    /// Applies operation op to the first qubit in the register
    /// # Input
    /// ## op
    /// An operation to be applied to the first qubit
    /// ## register 
    /// Qubit array to the first qubit of which the operation is applied
    operation ApplyToFirstQubit( op : (Qubit => ()), register : Qubit[] ) : () {
        body {
            op(register[0]);
        }
    }

    /// # See Also
    /// @"Microsoft.Quantum.Canon.ApplyToFirstQubit"
    operation ApplyToFirstQubitA( op : (Qubit => () : Adjoint), register : Qubit[] ) : () {
        body {
            op(register[0]);
        }
        adjoint auto
    }

    /// # See Also
    /// @"Microsoft.Quantum.Canon.ApplyToFirstQubit"
    operation ApplyToFirstQubitC( op : (Qubit => () : Controlled), register : Qubit[] ) : () {
        body {
            op(register[0]);
        }
        controlled auto
    }

    /// # See Also
    /// @"Microsoft.Quantum.Canon.ApplyToFirstQubit"
    operation ApplyToFirstQubitCA( op : (Qubit => () : Adjoint, Controlled), register : Qubit[] ) : () {
        body {
            op(register[0]);
        }
        adjoint auto
        controlled auto
        controlled adjoint auto
    }

    /// # Summary
    /// Applies operation `op` to the first two qubits in the register.
    ///
    /// # Input
    /// ## op
    /// An operation to be applied to the first two qubits
    /// ## register 
    /// Qubit array to the first two qubits of which the operation is applied.
    ///
    /// # Remarks
    /// This is equivalent to:
    /// ```qsharp
    /// op(register[0], register[1]);
    /// ```
    operation ApplyToFirstTwoQubits( op : ((Qubit,Qubit) => ()), register : Qubit[] ) : () {
        body {
            op(register[0],register[1]);
        }
    }

    /// # See Also
    /// @"Microsoft.Quantum.Canon.ApplyToFirstTwoQubits"
    operation ApplyToFirstTwoQubitsA( op : ((Qubit,Qubit) => () : Adjoint ), register : Qubit[] ) : () {
        body {
            op(register[0],register[1]);
        }
        adjoint auto
    }

    /// # See Also
    /// @"Microsoft.Quantum.Canon.ApplyToFirstTwoQubits"
    operation ApplyToFirstTwoQubitsC( op : ((Qubit,Qubit) => () : Controlled ), register : Qubit[] ) : () {
        body {
            op(register[0],register[1]);
        }
        controlled auto
    }

    /// # See Also
    /// @"Microsoft.Quantum.Canon.ApplyToFirstTwoQubits"
    operation ApplyToFirstTwoQubitsCA( op : ((Qubit,Qubit) => () : Adjoint, Controlled ), register : Qubit[] ) : () {
        body {
            op(register[0],register[1]);
        }
        adjoint auto
        controlled auto
        controlled adjoint auto
    }

    /// # Summary
    /// Applies operation `op` to the first three qubits in the register.
    /// # Input
    /// ## op
    /// An operation to be applied to the first three qubits
    /// ## register 
    /// Qubit array to the first three qubits of which the operation is applied.
    ///
    /// # Remarks
    /// This is equivalent to:
    /// ```qsharp
    /// op(register[0], register[1], register[2]);
    /// ```
    operation ApplyToFirstThreeQubits( op : ((Qubit,Qubit,Qubit) => ()), register : Qubit[] ) : () {
        body {
            op(register[0],register[1],register[2]);
        }
    }

    /// # See Also
    /// @"Microsoft.Quantum.Canon.ApplyToFirstThreeQubits"
    operation ApplyToFirstThreeQubitsA( op : ((Qubit,Qubit,Qubit) => () : Adjoint), register : Qubit[] ) : () {
        body {
            op(register[0],register[1],register[2]);
        }
        adjoint auto
    }

    /// # See Also
    /// @"Microsoft.Quantum.Canon.ApplyToFirstThreeQubits"
    operation ApplyToFirstThreeQubitsC( op : ((Qubit,Qubit,Qubit) => () : Controlled), register : Qubit[] ) : () {
        body {
            op(register[0],register[1],register[2]);
        }
        controlled auto
    }

    /// # See Also
    /// @"Microsoft.Quantum.Canon.ApplyToFirstThreeQubits"
    operation ApplyToFirstThreeQubitsCA( op : ((Qubit,Qubit,Qubit) => () : Adjoint, Controlled), register : Qubit[] ) : () {
        body {
            op(register[0],register[1],register[2]);
        }
        adjoint auto
        controlled auto
        controlled adjoint auto
    }
}