// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Diagnostics {
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Math;
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Logical;

    /// # Summary
    /// Represents the claim that a classical floating point value has the
    /// expected value up to a given
    /// absolute tolerance.
    ///
    /// # Input
    /// ## actual
    /// The number to be checked.
    /// ## expected
    /// The expected value.
    /// ## tolerance
    /// Absolute tolerance on the difference between actual and expected.
    function EqualityWithinToleranceFact(actual : Double, expected : Double, tolerance : Double) : Unit {
        let delta = actual - expected;
        if (delta > tolerance or delta < -tolerance) {
            FormattedFailure(actual, expected, "Values were not equal within tolerance.");
        }
    }

    /// # Summary
    /// Asserts that a classical floating point value has the expected value up to a
    /// small tolerance of 1e-10.
    ///
    /// # Input
    /// ## actual
    /// The number to be checked.
    /// ## expected
    /// The expected value.
    ///
    /// # Remarks
    /// This is equivalent to <xref:microsoft.quantum.diagnostics.equalitywithintolerancefact> with
    /// hardcoded tolerance of $10^{-10}$.
    function NearEqualityFactD(actual : Double, expected : Double) : Unit {
        EqualityWithinToleranceFact(actual, expected, 1e-10);
    }

    /// # Summary
    /// Asserts that a classical complex number has the expected value up to a
    /// small tolerance of 1e-10.
    ///
    /// # Input
    /// ## actual
    /// The number to be checked.
    /// ## expected
    /// The expected value.
    function NearEqualityFactC(actual : Complex, expected : Complex) : Unit {
        // Don't reduce to the base case of Fact, since we need to check two
        // conditions.
        let ((reA, imA), (reE, imE)) = (actual!, expected!);
        if (AbsD(reA - reE) >= 1e-10 or AbsD(imA - imE) >= 1e-10) {
            FormattedFailure(actual, expected, "Values were not equal within tolerance.");
        }
    }

    /// # Summary
    /// Asserts that a classical complex number has the expected value up to a
    /// small tolerance of 1e-10.
    ///
    /// # Input
    /// ## actual
    /// The number to be checked.
    /// ## expected
    /// The expected value.
    function NearEqualityFactCP(actual : ComplexPolar, expected : ComplexPolar) : Unit {
        return NearEqualityFactC(
            ComplexPolarAsComplex(actual),
            ComplexPolarAsComplex(expected)
        );
    }

    /// # Summary
    /// Asserts that a classical Int variable has the expected value.
    ///
    /// # Input
    /// ## actual
    /// The number to be checked.
    /// ## expected
    /// The expected value.
    ///
    /// ## message
    /// Failure message string to be used when the assertion is triggered.
    function EqualityFactI(actual : Int, expected : Int, message : String) : Unit {
        if (actual != expected) {
            FormattedFailure(actual, expected, message);
        }
    }

    /// # Summary
    /// Asserts that a classical BigInt variable has the expected value.
    ///
    /// # Input
    /// ## actual
    /// The number to be checked.
    /// ## expected
    /// The expected value.
    ///
    /// ## message
    /// Failure message string to be used when the assertion is triggered.
    function EqualityFactL(actual : BigInt, expected : BigInt, message : String) : Unit {
        if (actual != expected) {
            FormattedFailure(actual, expected, message);
        }
    }

    /// # Summary
    /// Asserts that a classical Bool variable has the expected value.
    ///
    /// # Input
    /// ## actual
    /// The variable to be checked.
    ///
    /// ## expected
    /// The expected value.
    ///
    /// ## message
    /// Failure message string to be used when the assertion is triggered.
    function EqualityFactB(actual : Bool, expected : Bool, message : String) : Unit {
        if (actual != expected) {
            FormattedFailure(actual, expected, message);
        }
    }

    /// # Summary
    /// Asserts that a classical Result variable has the expected value.
    ///
    /// # Input
    /// ## actual
    /// The variable to be checked.
    ///
    /// ## expected
    /// The expected value.
    ///
    /// ## message
    /// Failure message string to be used when the assertion is triggered.
    function EqualityFactR (actual : Result, expected : Result, message : String) : Unit {
        if (actual != expected) {
            FormattedFailure(actual, expected, message);
        }
    }

    /// # Summary
    /// Asserts that a complex number has the expected value.
    ///
    /// # Input
    /// ## actual
    /// The value to be checked.
    ///
    /// ## expected
    /// The expected value.
    ///
    /// ## message
    /// Failure message string to be used when the assertion is triggered.
    function EqualityFactC(actual : Complex, expected : Complex, message : String) : Unit {
        if (actual::Real != expected::Real or actual::Imag != expected::Imag) {
            FormattedFailure(actual, expected, message);
        }
    }

    /// # Summary
    /// Asserts that a complex number has the expected value.
    ///
    /// # Input
    /// ## actual
    /// The value to be checked.
    ///
    /// ## expected
    /// The expected value.
    ///
    /// ## message
    /// Failure message string to be used when the assertion is triggered.
    function EqualityFactCP(actual : ComplexPolar, expected : ComplexPolar, message : String) : Unit {
        let actualCartesian = ComplexPolarAsComplex(actual);
        let expectedCartesian = ComplexPolarAsComplex(expected);
        if (
            actualCartesian::Real != expectedCartesian::Real or
            actualCartesian::Imag != expectedCartesian::Imag
        ) {
            FormattedFailure(actual, expected, message);
        }
    }

    /// # Summary
    /// Asserts that two arrays of boolean values are equal.
    ///
    /// # Input
    /// ## actual
    /// The array that is produced by a test case of interest.
    /// ## expected
    /// The array that is expected from a test case of interest.
    /// ## message
    /// A message to be printed if the arrays are not equal.
    ///
    /// # See Also
    /// Microsoft.Quantum.Diagnostics.AllEqualityFactI
    function AllEqualityFactB(actual : Bool[], expected : Bool[], message : String) : Unit {
        let n = Length(actual);
        if (n != Length(expected)) {
            FormattedFailure(actual, expected, message);
        }

        Ignore(Mapped(EqualityFactB(_, _, message), Zip(actual, expected)));
    }

    /// # Summary
    /// Asserts that two arrays of integer values are equal.
    ///
    /// # Input
    /// ## actual
    /// The array that is produced by a test case of interest.
    /// ## expected
    /// The array that is expected from a test case of interest.
    /// ## message
    /// A message to be printed if the arrays are not equal.
    ///
    /// # See Also
    /// Microsoft.Quantum.Diagnostics.AllEqualityFactB
    function AllEqualityFactI(actual : Int[], expected : Int[], message : String) : Unit {
        let n = Length(actual);
        if (n != Length(expected)) {
            FormattedFailure(actual, expected, message);
        }

        Ignore(Mapped(EqualityFactI(_, _, message), Zip(actual, expected)));
    }

}
