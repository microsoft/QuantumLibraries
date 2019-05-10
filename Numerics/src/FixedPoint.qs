// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Arithmetic {
    open Microsoft.Quantum.Convert;
    open Microsoft.Quantum.Math;
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Measurement;
    open Microsoft.Quantum.Primitive;
	open Microsoft.Quantum.Diagnostics;

    /// # Summary
    /// FixedPoint type. Consists of an integer that is equal to the number of
    /// qubits to the left of the binary point, i.e., qubits of weight greater
    /// than or equal to 1, and a quantum register.
    newtype FixedPoint = (Int, Qubit[]);

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

    /// # Summary
    /// Initialize a quantum fixed-point number to a classical constant.
    ///
    /// # Input
    /// ## constant
    /// Constant to which to initialize the quantum fixed-point number.
    /// ## fp
    /// Fixed-point number (of type FixedPoint) to initialize.
    operation FixedPointInit(constant : Double, fp : FixedPoint) : Unit{
        body (...) {
            let (p, q) = fp!;
            let n = Length(q);
            let sign = constant < 0.;
            mutable rescaledConstant = PowD(2., IntAsDouble(n-p)) * AbsD(constant) + 0.5;
            mutable keepAdding = sign;
            for (i in 0..n-1) {
                let intConstant = Floor(rescaledConstant);
                set rescaledConstant = 0.5 * rescaledConstant;
                mutable currentBit = (intConstant &&& 1) == (sign ? 0 | 1);
                if (keepAdding) {
                    set keepAdding = currentBit;
                    set currentBit = not currentBit;
                }
                if (currentBit) {
                    X(q[i]);
                }
            }
        }
        controlled auto;
        adjoint self;
        adjoint controlled auto;
    }

    /// # Summary
    /// Measure a fixed-point number and return its value as Double.
    ///
    /// # Input
    /// ## fp
    /// Fixed-point number to measure.
    operation MeasureFixedPoint(fp : FixedPoint) : Double {
        let (p, xs) = fp!;
        let n = Length(xs);
        let sign = MResetZ(xs[n-1]) == One;
        mutable keepAdding = sign;
        mutable fpAsDouble = 0.;
        for (i in 0..n-2) {
            mutable currentRes = MResetZ(xs[i]) == (sign ? Zero | One);
            if (keepAdding) {
                set keepAdding = currentRes;
                set currentRes = not currentRes;
            }
            set fpAsDouble = fpAsDouble * 0.5 + (currentRes == true ? 1. | 0.);
        }
        return (sign ? -1.0 | 1.0) * fpAsDouble * PowD(2.0, IntAsDouble(p-2));
    }

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
                WithCA(FixedPointInit(constant, _), FixedPointAddition(_, fp), tmpFp);
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
    operation FixedPointGreaterThan(fp1 : FixedPoint, fp2 : FixedPoint,
                                    result : Qubit) : Unit {
        body(...) {
            let (px, xs) = fp1!;
            let (py, ys) = fp2!;

            IdenticalFormatFactFP([fp1, fp2]);
            SignedIntegerGreaterThan(SignedLittleEndian(LittleEndian(xs)),
                                     SignedLittleEndian(LittleEndian(ys)),
                                     result);
        }
        adjoint auto;
        controlled auto;
        adjoint controlled auto;
    }

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
    operation FixedPointMultiplication(fp1 : FixedPoint, fp2 : FixedPoint,
                                       result : FixedPoint) : Unit {
        body(...) {
            (Controlled FixedPointMultiplication) (new Qubit[0],
                                                   (fp1, fp2, result));
        }
        controlled (controls, ...){
            IdenticalFormatFactFP([fp1, fp2, result]);
            AssertAllZeroFP(result);
            let (px, xs) = fp1!;
            let (py, ys) = fp2!;
            let (pz, zs) = result!;
            let n = Length(xs);

            using (tmpResult = Qubit[2*n]){
                let xsInt = SignedLittleEndian(LittleEndian(xs));
                let ysInt = SignedLittleEndian(LittleEndian(ys));
                let tmpResultInt = SignedLittleEndian(
                    LittleEndian(tmpResult));
                SignedIntegerMultiplication(xsInt, ysInt,
                                            tmpResultInt);
                (Controlled ApplyToEachCA)(controls,
                                           (CNOT,
                                            Zip(tmpResult[n-px..2*n-px-1], zs)));
                (Adjoint SignedIntegerMultiplication)(xsInt, ysInt,
                                                      tmpResultInt);
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
    operation FixedPointSquare(fp : FixedPoint, result : FixedPoint) : Unit {
        body(...) {
            (Controlled FixedPointSquare) (new Qubit[0],
                                           (fp, result));
        }
        controlled (controls, ...){
            IdenticalFormatFactFP([fp, result]);
            AssertAllZeroFP(result);
            let (px, xs) = fp!;
            let (py, ys) = result!;
            let n = Length(xs);

            using (tmpResult = Qubit[2*n]){
                let xsInt = SignedLittleEndian(LittleEndian(xs));
                let tmpResultInt = SignedLittleEndian(
                    LittleEndian(tmpResult));
                SignedIntegerSquare(xsInt, tmpResultInt);
                (Controlled ApplyToEachCA)(controls,
                                           (CNOT,
                                            Zip(tmpResult[n-px..2*n-px-1], ys)));
                (Adjoint SignedIntegerSquare)(xsInt, tmpResultInt);
            }
        }
        adjoint auto;
        adjoint controlled auto;
    }
    
    /// # Summary
    /// Polynomial evaluation in a fixed-point representation
    ///
    /// # Input
    /// ## coefficients
    /// Coefficients of the polynomial as a double array, i.e., the array
    /// $[a_0, a_1, ..., a_d]$ for the polynomial
    /// $P(x) = a_0 + a_1 x + \cdots + a_d x^d$.
    /// ## fpx
    /// Input FixedPoint number for which to evaluate the polynomial.
    /// ## result
    /// Output FixedPoint number which will hold P(x). Must be in state
    /// $\ket{0}$ initially.
    operation FixedPointPolynomial(coefficients : Double[], fpx : FixedPoint,
                                   result : FixedPoint) : Unit {
        body (...) {
            (Controlled FixedPointPolynomial) (new Qubit[0],
                (coefficients, fpx, result));
        }
        controlled (controls, ...) {
            IdenticalFormatFactFP([fpx, result]);
            AssertAllZeroFP(result);
            let degree = Length(coefficients) - 1;
            let (p, q) = fpx!;
            let n = Length(q);
            if (degree == 0){
                (Controlled FixedPointInit)(controls,
                    (coefficients[0], result));
            }
            elif (degree > 0) {
                // initialize ancillary register to a_d
                using (qubits = Qubit[n * degree]){
                    let firstIterate = FixedPoint(p,
                        qubits[(degree-1)*n..degree*n-1]);
                    (Controlled FixedPointInit)(controls,
                        (coefficients[degree], firstIterate));
                    for (d in degree..(-1)..2) {
                        let currentIterate = FixedPoint(p, qubits[(d-1)*n..d*n-1]);
                        let nextIterate = FixedPoint(p, qubits[(d-2)*n..(d-1)*n-1]);
                        // multiply by x and then add current coefficient
                        FixedPointMultiplication(currentIterate, fpx, nextIterate);
                        (Controlled FixedPointAdditionConstant)(controls,
                            (coefficients[d-1], nextIterate));
                    }
                    let finalIterate = FixedPoint(p, qubits[0..n-1]);
                    // final multiplication into the result register
                    FixedPointMultiplication(finalIterate, fpx, result);
                    // add a_0 to complete polynomial evaluation and
                    (Controlled FixedPointAdditionConstant)(controls,
                        (coefficients[0], result));
                    // uncompute intermediate results
                    for (d in 2..degree) {
                        let currentIterate = FixedPoint(p, qubits[(d-1)*n..d*n-1]);
                        let nextIterate = FixedPoint(p, qubits[(d-2)*n..(d-1)*n-1]);
                        (Controlled Adjoint FixedPointAdditionConstant)(controls,
                            (coefficients[d-1], nextIterate));
                        (Adjoint FixedPointMultiplication)(currentIterate, fpx,
                                                           nextIterate);
                    }
                    (Controlled FixedPointInit)(controls,
                        (coefficients[degree], firstIterate));
                }
            }
        }
        adjoint auto;
        adjoint controlled auto;
    }

    /// # Summary
    /// Evaluation of an even polynomial in a fixed-point representation
    ///
    /// # Input
    /// ## coefficients
    /// Coefficients of the even polynomial as a double array, i.e., the array
    /// $[a_0, a_1, ..., a_k]$ for the even polynomial
    /// $P(x) = a_0 + a_1 x^2 + \cdots + a_k x^{2k}$.
    /// ## fpx
    /// Input FixedPoint number for which to evaluate the polynomial.
    /// ## result
    /// Output FixedPoint number which will hold P(x). Must be in state
    /// $\ket{0}$ initially.
    operation FixedPointEvenPolynomial(coefficients : Double[], fpx : FixedPoint,
                                       result : FixedPoint) : Unit {
        body (...) {
            (Controlled FixedPointEvenPolynomial) (new Qubit[0],
                (coefficients, fpx, result));
        }
        controlled (controls, ...) {
            IdenticalFormatFactFP([fpx, result]);
            AssertAllZeroFP(result);
            let halfDegree = Length(coefficients) - 1;
            let (p, q) = fpx!;
            let n = Length(q);

            if (halfDegree == 0){
                (Controlled FixedPointInit)(controls,
                    (coefficients[0], result));
            }
            elif (halfDegree > 0) {
                // initialize ancillary register to a_d
                using (xsSquared = Qubit[n]){
                    let fpxSquared = FixedPoint(p, xsSquared);
					WithCA(FixedPointSquare(fpx, _),
						(Controlled FixedPointPolynomial)(controls,
							(coefficients, _, result)),
						fpxSquared);
                }
            }
        }
        adjoint auto;
        adjoint controlled auto;
    }

    /// # Summary
    /// Evaluation of an odd polynomial in a fixed-point representation
    ///
    /// # Input
    /// ## coefficients
    /// Coefficients of the odd polynomial as a double array, i.e., the array
    /// $[a_0, a_1, ..., a_k]$ for the odd polynomial
    /// $P(x) = a_0 x + a_1 x^3 + \cdots + a_k x^{2k+1}$.
    /// ## fpx
    /// Input FixedPoint number for which to evaluate the polynomial.
    /// ## result
    /// Output FixedPoint number which will hold P(x). Must be in state
    /// $\ket{0}$ initially.
    operation FixedPointOddPolynomial(coefficients : Double[], fpx : FixedPoint,
                                      result : FixedPoint) : Unit {
        body (...) {
            (Controlled FixedPointOddPolynomial) (new Qubit[0],
                (coefficients, fpx, result));
        }
        controlled (controls, ...) {
            IdenticalFormatFactFP([fpx, result]);
            AssertAllZeroFP(result);
            let halfDegree = Length(coefficients) - 1;
            let (p, q) = fpx!;
            let n = Length(q);
            if (halfDegree >= 0) {
                using (tmpResult = Qubit[n]) {
                    let tmpResultFp = FixedPoint(p, tmpResult);
					WithCA(FixedPointEvenPolynomial(coefficients, _, _),
						   (Controlled FixedPointMultiplication)(controls,
																 (_, _, result)),
						   (fpx, tmpResultFp));
                }
            }
        }
        adjoint auto;
        adjoint controlled auto;
    }

    /// # Summary
    /// Computes 1/x for a fixed-point number x.
    ///
    /// # Input
    /// ## x
    /// FixedPoint number to invert
    /// ## result
    /// FixedPoint number that will hold the result. Must be initialized to 0.
    operation FixedPointReciprocal(x : FixedPoint, result : FixedPoint) : Unit {
        body (...) {
            (Controlled FixedPointReciprocal) (new Qubit[0], (x, result));
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
                (Controlled IntegerInversion2s)
                    ([sign], SignedLittleEndian(LittleEndian(xs)));
                IntegerReciprocal(LittleEndian(xs), LittleEndian(tmpRes));
                (Controlled ApplyToEachCA)(controls,
                    (CNOT, Zip(tmpRes[p+pRes-1+n-Length(rs)..Min([n+p+pRes, 2*n-1])], rs)));
                (Controlled IntegerInversion2s)([sign], SignedLittleEndian(LittleEndian(rs)));
                (Adjoint IntegerReciprocal)(LittleEndian(xs), LittleEndian(tmpRes));
                (Controlled Adjoint IntegerInversion2s)
                    ([sign], SignedLittleEndian(LittleEndian(xs)));
                CNOT(Tail(xs), sign);
            }
        }
        adjoint auto;
        adjoint controlled auto;
    }
}