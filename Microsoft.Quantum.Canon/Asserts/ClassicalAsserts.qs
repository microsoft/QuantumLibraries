// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {

    /// # Summary
    /// Asserts that a classical floating point variable has the expected value up to a given
    /// absolute tolerance.
    ///
    /// # Input
    /// ## actual
    /// The number to be checked.
    ///
    /// ## expected
    /// The expected value.
    ///
    /// ## tolerance
    /// Absolute tolerance on the difference between actual and expected.
    function AssertAlmostEqualTol(actual : Double, expected : Double, tolerance : Double) : () {
        let delta = actual - expected;
        // TODO: rewrite using Abs.
        if (delta > tolerance || delta < -tolerance) {
            // TODO: use string interpolation here.
            fail $"Assertion failed.\n\tExpected: {expected}.\n\tActual:   {actual}";
        }
    }

    /// # Summary
    /// Asserts that a classical floating point variable has the expected value up to a
    /// small tolerance of 1e-10.
    ///
    /// # Input
    /// ## actual
    /// The number to be checked.
    ///
    /// ## expected
    /// The expected value.
    ///
    /// # Remarks
    /// This is equivalent to @"Microsoft.Quantum.Canon.AssertAlmostEqualTol" with
    /// hardcoded tolerance=1e-10.
    function AssertAlmostEqual(actual : Double, expected : Double) : () {
        AssertAlmostEqualTol(actual, expected, 1e-10);
    }

    // FIXME: the following asserts should be made generic if possible.

    /// # Summary
    /// Asserts that a classical Int variable has the expected value.
    ///
    /// # Input
    /// ## actual
    /// The number to be checked.
    ///
    /// ## expected
    /// The expected value.
    ///
    /// ## message
    /// Failure message string to be used when the assertion is triggered.
    function AssertIntEqual ( actual : Int, expected : Int, message : String ) : () {
        if ( actual != expected )
        {
            fail message;
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
    function AssertBoolEqual ( actual : Bool, expected : Bool, message : String ) : () {
        if( actual != expected ) {
            fail message;
        }
    }

	function AssertBoolArrayEqual ( actual : Bool[], expected : Bool[], message : String ) : () {
		let n = Length(actual); 
		if (n != Length(expected)) {
			fail message;
		}	
		for (idx in 0..(n-1)) {
			if( actual[idx] != expected[idx] ) {
				fail message;
			}
        }
    }

}
