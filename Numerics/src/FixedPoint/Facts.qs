// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Arithmetic {
	open Microsoft.Quantum.Diagnostics;

    /// # Summary
    /// Helper function to assert that a quantum fixed-point number is
    /// initialized to zero, i.e., all qubits are in state $\ket{0}$.
    operation AssertAllZeroFP(fp : FixedPoint) : Unit {
        body (...) {
            let (p, xs) = fp!;
            AssertAllZero(xs);
        }
        adjoint auto;
        controlled auto;
        adjoint controlled auto;
    }

    /// # Summary
    /// Assert that all fixed-point numbers in the provided array
    /// have identical point positions and qubit numbers.
    ///
    /// # Input
    /// ## fixedPoints
    /// Array of quantum fixed-point numbers that will be checked for
    /// compatibility (using assertions).
    function IdenticalFormatFactFP(fixedPoints : FixedPoint[]) : Unit {
        if (Length(fixedPoints) == 0){
            return ();
        }
        let (position, register) = fixedPoints[0]!;
        EqualityFactB(position > 0, true,
            "Point position must be greater than zero.");
        let n = Length(register);
        for (fp in fixedPoints[1..Length(fixedPoints)-1]) {
            let (pos, reg) = fp!;
            EqualityFactI(pos, position,
                "FixedPoint numbers must have identical binary point position.");
            EqualityFactI(Length(reg), n,
                "FixedPoint numbers must have identical number of qubits.");
        }
    }

    /// # Summary
    /// Assert that all fixed-point numbers in the provided array
    /// have identical point positions when counting from the least-
    /// significant bit. I.e., number of bits minus point position must
    /// be constant for all fixed-point numbers in the array.
    ///
    /// # Input
    /// ## fixedPoints
    /// Array of quantum fixed-point numbers that will be checked for
    /// compatibility (using assertions).
    function IdenticalPointPosFactFP(fixedPoints : FixedPoint[]) : Unit {
        if (Length(fixedPoints) == 0){
            return ();
        }
        let (position, register) = fixedPoints[0]!;
        EqualityFactB(position > 0, true,
            "Point position must be greater than zero.");
        let n = Length(register);
        for (fp in fixedPoints[1..Length(fixedPoints)-1]) {
            let (pos, reg) = fp!;
            EqualityFactI(Length(reg)-pos, n-position,
                "FixedPoint numbers must have identical point alignment.");
        }
    }
}