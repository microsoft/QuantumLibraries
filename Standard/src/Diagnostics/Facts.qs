// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Diagnostics {
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Arrays;

    /// # Summary
    /// Declares that a classical condition is true.
    ///
    /// # Input
    /// ## actual
    /// The condition to be declared.
    /// ## message
    /// Failure message string to be printed in the case that the classical
    /// condition is false.
    function Fact(actual : Bool, message : String) : Unit {
        if (not actual) { fail message; }
    }

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
            fail $"Fact was false. Expected: '{expected}'. Actual: '{actual}'";
        }
    }

    /// # Summary
    /// Asserts that a classical floating point variable has the expected value up to a
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
    function NearEqualityFact(actual : Double, expected : Double) : Unit {
        EqualityWithinToleranceFact(actual, expected, 1E-10);
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
            fail $"{actual} ≠ {expected}: {message}";
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
            fail $"{actual} ≠ {expected}: {message}";
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
            fail $"{actual} ≠ {expected}: {message}";
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
            fail $"{actual} ≠ {expected}: {message}";
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
    function AllEqualityFactB(actual : Bool[], expected : Bool[], message : String) : Unit {
        let n = Length(actual);
        if (n != Length(expected)) {
            fail message;
        }

        Ignore(Mapped(EqualityFactB(_, _, message), Zip(actual, expected)));
    }

}
