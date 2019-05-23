// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Arithmetic {
    open Microsoft.Quantum.Canon;

    /// # Summary
    /// Adds a classical constant to a quantum fixed-point number.
    ///
    /// # Input
    /// ## constant
    /// Constant to add to the quantum fixed-point number.
    /// ## fp
    /// Fixed-point number to which the constant will
    /// be added.
    operation AddConstantFxP(constant : Double, fp : FixedPoint) : Unit is Adj + Ctl {
        let (px, xs) = fp!;
        let n = Length(xs);
        using (ys = Qubit[n]){
            let tmpFp = FixedPoint(px, ys);
            ApplyWithCA(InitFxP(constant, _), AddFxP(_, fp), tmpFp);
        }
    }

    /// # Summary
    /// Adds two fixed-point numbers stored in quantum registers.
    ///
    /// # Input
    /// ## fp1
    /// First fixed-point number
    /// ## fp2
    /// Second fixed-point number, will be updated to contain the sum of the
    /// two inputs.
    ///
    /// # Remarks
    /// The current implementation requires the two fixed-point numbers
    /// to have the same point position counting from the least-significant
    /// bit, i.e., $n_i$ and $p_i$ must be equal.
    operation AddFxP(fp1 : FixedPoint, fp2 : FixedPoint) : Unit is Adj + Ctl {
        let (px, xs) = fp1!;
        let (py, ys) = fp2!;

        IdenticalPointPosFactFxP([fp1, fp2]);

        AddI(LittleEndian(xs), LittleEndian(ys));
    }
}
