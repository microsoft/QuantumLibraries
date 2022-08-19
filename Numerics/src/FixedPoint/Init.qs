// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Arithmetic {
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Convert;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Math;

    /// # Summary
    /// Initialize a quantum fixed-point number to a classical constant.
    ///
    /// # Input
    /// ## constant
    /// Constant to which to initialize the quantum fixed-point number.
    /// ## fp
    /// Fixed-point number (of type FixedPoint) to initialize.
    operation PrepareFxP(constant : Double, fp : FixedPoint)
    : Unit is Adj + Ctl {
        let bits = FixedPointAsBoolArray(fp::IntegerBits, Length(fp::Register) - fp::IntegerBits, constant);
        ApplyPauliFromBitString(PauliX, true, bits, fp::Register);
    }
}
