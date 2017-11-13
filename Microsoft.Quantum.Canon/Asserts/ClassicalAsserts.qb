// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {

	function AssertAlmostEqualTol(actual : Double, expected : Double, tolerance : Double) : () {
		let delta = actual - expected;
		// TODO: rewrite using Abs.
		if (delta > tolerance || delta < -tolerance) {
			// TODO: use string interpolation here.
			fail "Assertion failed.";
		}
	}

	function AssertAlmostEqual(actual : Double, expected : Double) : () {
		AssertAlmostEqualTol(actual, expected, 1e-10);
	}

    // FIXME: the following asserts should be made generic if possible.

    function AssertIntEqual ( actual : Int, expected : Int, message : String ) : () {
        if ( actual != expected )
        {
            fail message;
        }
    }

	function AssertBoolEqual ( actual : Bool, expected : Bool, message : String ) : () {
        if( actual != expected ) {
            fail message;
        }
    }

}
