// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Arithmetic {

    /// # Summary
    /// Comparison of two fixed-point numbers
    ///
    /// # Input
    /// ## fp1
    /// First fixed-point number (of type FixedPoint)
    /// ## fp2
    /// Second fixed-point number (of type FixedPoint)
    /// ## result
    /// Result of the comparison. Will be flipped if `fp1 > fp2`.
    ///
    /// # Remarks
    /// The current implementation requires the two fixed-point numbers
    /// to have the same point position and the same number of qubits.
    operation CompareGTFxP(fp1 : FixedPoint, fp2 : FixedPoint,
                                    result : Qubit) : Unit {
        body(...) {
            let (px, xs) = fp1!;
            let (py, ys) = fp2!;

            IdenticalFormatFactFxP([fp1, fp2]);
            CompareGTSI(SignedLittleEndian(LittleEndian(xs)),
                        SignedLittleEndian(LittleEndian(ys)),
                        result);
        }
        adjoint auto;
        controlled auto;
        adjoint controlled auto;
    }
}