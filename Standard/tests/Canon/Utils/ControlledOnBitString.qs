// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Tests {
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Diagnostics;
    open Microsoft.Quantum.Arrays;

    /// # Summary
    /// Applies an Ry rotation controlled on
    /// the state $\ket{001}$, leaving additional
    /// control qubits uncontrolled.
    ///
    /// # Input
    /// ## register
    /// A register of at least four qubits, the last of which is considered to
    /// be the target.
    internal operation ApplyRyControlledOn001(register : Qubit[])
    : Unit is Adj + Ctl {
        let controls = Most(register);
        let target = Tail(register);

        within {
            X(controls[0]);
            X(controls[1]);
        } apply {
            Controlled Ry(controls[...2], (-1.234, target));
        }
    }

    internal operation ApplyControlledOnBitStringToFlatRegister(
        op : ((Qubit[], Qubit) => Unit is Adj + Ctl),
        register : Qubit[]
    ) : Unit is Adj + Ctl {
        op(Most(register), Tail(register));
    }

    @Test("QuantumSimulator")
    operation CheckControlledOnBitString() : Unit {
        let controlledOp = ControlledOnBitString(
            [false, false, true],
            Ry(-1.234, _)
        );
        AssertOperationsEqualReferenced(4,
            ApplyRyControlledOn001,
            controlledOp
        );

        // Check that extra qubits are uncontrolled.
        AssertOperationsEqualReferenced(5,
            ApplyRyControlledOn001,
            controlledOp
        );
    }

}
