// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Arithmetic {
    open Microsoft.Quantum.Canon;
    
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
					ApplyWithCA(FixedPointSquare(fpx, _),
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
					ApplyWithCA(FixedPointEvenPolynomial(coefficients, _, _),
						   (Controlled FixedPointMultiplication)(controls,
																 (_, _, result)),
						   (fpx, tmpResultFp));
                }
            }
        }
        adjoint auto;
        adjoint controlled auto;
    }
}