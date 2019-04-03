// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {
    open Microsoft.Quantum.Primitive;
    open Microsoft.Quantum.Diagnostics;

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Diagnostics.ClaimEqualWithinTolerance" instead.
    function AssertAlmostEqualTol (actual : Double, expected : Double, tolerance : Double) : Unit {
        Renamed("Microsoft.Quantum.Canon.AssertAlmostEqualTol", "Microsoft.Quantum.Diagnostics.ClaimEqualWithinTolerance");
        return ClaimEqualWithinTolerance(actual, expected, tolerance);
    }

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Diagnostics.ClaimAlmostEqual" instead.
    function AssertAlmostEqual(actual : Double, expected : Double) : Unit {
        Renamed("Microsoft.Quantum.Canon.AssertAlmostEqual", "Microsoft.Quantum.Diagnostics.ClaimAlmostEqual");
        AssertAlmostEqual(actual, expected);
    }

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Diagnostics.ClaimEqualI" instead.
    function AssertIntEqual (actual : Int, expected : Int, message : String) : Unit {
        Renamed("Microsoft.Quantum.Canon.AssertIntEqual", "Microsoft.Quantum.Diagnostics.ClaimEqualI");
        ClaimEqualI(actual, expected, message);
    }

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Diagnostics.ClaimEqualB" instead.
    function AssertBoolEqual(actual : Bool, expected : Bool, message : String) : Unit {
        Renamed("Microsoft.Quantum.Canon.AssertBoolEqual", "Microsoft.Quantum.Diagnostics.ClaimEqualB");
        ClaimEqualB(actual, expected, message);
    }

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Diagnostics.ClaimEqualR" instead.
    function AssertResultEqual(actual : Result, expected : Result, message : String) : Unit {
        Renamed("Microsoft.Quantum.Canon.AssertResultEqual", "Microsoft.Quantum.Diagnostics.ClaimEqualR");
        ClaimEqualR(actual, expected, message);
    }

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Diagnostics.ClaimAllEqualB" instead.
    function AssertBoolArrayEqual (actual : Bool[], expected : Bool[], message : String) : Unit {
        Renamed("Microsoft.Quantum.Canon.AssertBoolArrayEqual", "Microsoft.Quantum.Diagnostics.ClaimAllEqualB");
        ClaimAllEqualB(actual, expected, message);
    }

}
