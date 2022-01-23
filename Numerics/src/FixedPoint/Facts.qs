// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Arithmetic {
    open Microsoft.Quantum.Diagnostics;
    open Microsoft.Quantum.Arrays;

    /// # Summary
    /// Asserts that a quantum fixed-point number is
    /// initialized to zero.
    ///
    /// # Description
    /// This assertion succeeds when all qubits are in state $\ket{0}$,
    /// representing that the register encodes the fixed-point number $0.0$.
    operation AssertAllZeroFxP(fp : FixedPoint) : Unit is Adj + Ctl {
        let (p, xs) = fp!;
        AssertAllZero(xs);
    }

    /// # Summary
    /// Assert that all fixed-point numbers in the provided array
    /// have identical point positions and qubit numbers.
    ///
    /// # Input
    /// ## fixedPoints
    /// Array of quantum fixed-point numbers that will be checked for
    /// compatibility (using assertions).
    function IdenticalFormatFactFxP(fixedPoints : FixedPoint[]) : Unit {
        if (Length(fixedPoints) == 0) {
            return ();
        }
        let (position, register) = fixedPoints[0]!;
        Fact(position > 0, "Point position must be greater than zero.");
        let n = Length(register);
        for fp in Most(fixedPoints) {
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
    function IdenticalPointPosFactFxP(fixedPoints : FixedPoint[]) : Unit {
        if (Length(fixedPoints) == 0) {
            return ();
        }
        let (position, register) = fixedPoints[0]!;
        Fact(position > 0, "Point position must be greater than zero.");
        let n = Length(register);
        for fp in Most(fixedPoints) {
            let (pos, reg) = fp!;
            EqualityFactI(Length(reg) - pos, n - position,
                "FixedPoint numbers must have identical point alignment.");
        }
    }
}
