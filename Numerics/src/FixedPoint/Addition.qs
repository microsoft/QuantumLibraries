// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Arithmetic {
    open Microsoft.Quantum.Canon;

    /// # Summary
    /// Addition of a classical constant to a quantum fixed-point number
    ///
    /// # Input
    /// ## constant
    /// Constant to add to the quantum fixed-point number.
    /// ## fp
    /// Fixed-point number (of type FixedPoint), to which the constant will
    /// be added.
    operation FixedPointAdditionConstant(constant : Double, fp : FixedPoint) : Unit {
        body(...) {
            let (px, xs) = fp!;
            let n = Length(xs);
            using (ys = Qubit[n]){
                let tmpFp = FixedPoint(px, ys);
                ApplyWithCA(FixedPointInit(constant, _), FixedPointAddition(_, fp), tmpFp);
            }
        }
        adjoint auto;
        controlled auto;
        adjoint controlled auto;
    }

    /// # Summary
    /// Addition of two fixed-point numbers
    ///
    /// # Input
    /// ## fp1
    /// First fixed-point number (of type FixedPoint)
    /// ## fp2
    /// Second fixed-point number (of type FixedPoint), will be updated
    /// to contain the sum of the two inputs.
    ///
    /// # Remarks
    /// The current implementation requires the two fixed-point numbers
    /// to have the same point position counting from the least-significant
    /// bit, i.e., n_i - p_i must be equal.
    operation FixedPointAddition(fp1 : FixedPoint, fp2 : FixedPoint) : Unit {
        body(...) {
            let (px, xs) = fp1!;
            let (py, ys) = fp2!;

            IdenticalPointPosFactFP([fp1, fp2]);

            IntegerAddition(LittleEndian(xs), LittleEndian(ys));
        }
        adjoint auto;
        controlled auto;
        adjoint controlled auto;
    }
}