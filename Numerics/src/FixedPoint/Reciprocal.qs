// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Arithmetic {
    open Microsoft.Quantum.Math;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Diagnostics;

    /// # Summary
    /// Computes 1/x for a fixed-point number x.
    ///
    /// # Input
    /// ## x
    /// FixedPoint number to invert
    /// ## result
    /// FixedPoint number that will hold the result. Must be initialized to 0.
    operation ComputeReciprocalFxP(x : FixedPoint, result : FixedPoint) : Unit {
        body (...) {
            (Controlled ComputeReciprocalFxP) (new Qubit[0], (x, result));
        }
        controlled (controls, ...) {
            let (p, xs) = x!;
            let (pRes, rs) = result!;
            let n = Length(xs);
            AssertAllZero(rs);
            EqualityFactB(p+pRes-1+n >= Length(rs), true,
                            "Output register is too wide.");
            using ((sign, tmpRes) = (Qubit(), Qubit[2*n])) {
                CNOT(Tail(xs), sign);
                (Controlled Invert2sI)
                    ([sign], SignedLittleEndian(LittleEndian(xs)));
                ComputeReciprocalI(LittleEndian(xs), LittleEndian(tmpRes));
                (Controlled ApplyToEachCA)(controls,
                    (CNOT, Zip(tmpRes[p+pRes-1+n-Length(rs)..Min([n+p+pRes, 2*n-1])], rs)));
                (Controlled Invert2sI)([sign], SignedLittleEndian(LittleEndian(rs)));
                (Adjoint ComputeReciprocalI)(LittleEndian(xs), LittleEndian(tmpRes));
                (Controlled Adjoint Invert2sI)
                    ([sign], SignedLittleEndian(LittleEndian(xs)));
                CNOT(Tail(xs), sign);
            }
        }
        adjoint auto;
        adjoint controlled auto;
    }
}