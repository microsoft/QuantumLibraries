// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Arithmetic {
	open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Canon;

    /// # Summary
    /// Multiplication of two fixed-point numbers
    ///
    /// # Input
    /// ## fp1
    /// First fixed-point number (of type FixedPoint)
    /// ## fp2
    /// Second fixed-point number (of type FixedPoint)
    /// ## result
    /// Result fixed-point number (of type FixedPoint),
    /// must be in state $\ket{0}$ initially.
    ///
    /// # Remarks
    /// The current implementation requires the three fixed-point numbers
    /// to have the same point position and the same number of qubits.
    operation MultiplyFxP(fp1 : FixedPoint, fp2 : FixedPoint,
                                       result : FixedPoint) : Unit {
        body(...) {
            (Controlled MultiplyFxP) (new Qubit[0],
                                                   (fp1, fp2, result));
        }
        controlled (controls, ...){
            IdenticalFormatFactFxP([fp1, fp2, result]);
            AssertAllZeroFxP(result);
            let (px, xs) = fp1!;
            let (py, ys) = fp2!;
            let (pz, zs) = result!;
            let n = Length(xs);

            using (tmpResult = Qubit[2*n]){
                let xsInt = SignedLittleEndian(LittleEndian(xs));
                let ysInt = SignedLittleEndian(LittleEndian(ys));
                let tmpResultInt = SignedLittleEndian(
                    LittleEndian(tmpResult));
                MultiplySI(xsInt, ysInt, tmpResultInt);
                (Controlled ApplyToEachCA)(controls,
                                           (CNOT,
                                            Zip(tmpResult[n-px..2*n-px-1], zs)));
                (Adjoint MultiplySI)(xsInt, ysInt, tmpResultInt);
            }
        }
        adjoint auto;
        adjoint controlled auto;
    }

    /// # Summary
    /// Square a fixed-point number
    ///
    /// # Input
    /// ## fp
    /// Fixed-point number (of type FixedPoint)
    /// ## result
    /// Result fixed-point number (of type FixedPoint),
    /// must be in state $\ket{0}$ initially.
    operation SquareFxP(fp : FixedPoint, result : FixedPoint) : Unit {
        body(...) {
            (Controlled SquareFxP) (new Qubit[0],
                                           (fp, result));
        }
        controlled (controls, ...){
            IdenticalFormatFactFxP([fp, result]);
            AssertAllZeroFxP(result);
            let (px, xs) = fp!;
            let (py, ys) = result!;
            let n = Length(xs);

            using (tmpResult = Qubit[2*n]){
                let xsInt = SignedLittleEndian(LittleEndian(xs));
                let tmpResultInt = SignedLittleEndian(
                    LittleEndian(tmpResult));
                SquareSI(xsInt, tmpResultInt);
                (Controlled ApplyToEachCA)(controls,
                                           (CNOT,
                                            Zip(tmpResult[n-px..2*n-px-1], ys)));
                (Adjoint SquareSI)(xsInt, tmpResultInt);
            }
        }
        adjoint auto;
        adjoint controlled auto;
    }
}