// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Diagnostics;

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Diagnostics.EqualityWithinToleranceFact" instead.
    @Deprecated("Microsoft.Quantum.Diagnostics.EqualityWithinToleranceFact")
    function AssertAlmostEqualTol (actual : Double, expected : Double, tolerance : Double) : Unit {
        return EqualityWithinToleranceFact(actual, expected, tolerance);
    }

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Diagnostics.NearEqualityFactD" instead.
    @Deprecated("Microsoft.Quantum.Diagnostics.NearEqualityFactD")
    function AssertAlmostEqual(actual : Double, expected : Double) : Unit {
        NearEqualityFactD(actual, expected);
    }

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Diagnostics.NearEqualityFactD" instead.
    @Deprecated("Microsoft.Quantum.Diagnostics.NearEqualityFactD")
    function NearEqualityFact(actual : Double, expected : Double) : Unit {
        NearEqualityFactD(actual, expected);
    }

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Diagnostics.EqualityFactI" instead.
    @Deprecated("Microsoft.Quantum.Diagnostics.EqualityFactI")
    function AssertIntEqual (actual : Int, expected : Int, message : String) : Unit {
        EqualityFactI(actual, expected, message);
    }

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Diagnostics.EqualityFactB" instead.
    @Deprecated("Microsoft.Quantum.Diagnostics.EqualityFactB")
    function AssertBoolEqual(actual : Bool, expected : Bool, message : String) : Unit {
        EqualityFactB(actual, expected, message);
    }

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Diagnostics.EqualityFactR" instead.
    @Deprecated("Microsoft.Quantum.Diagnostics.EqualityFactR")
    function AssertResultEqual(actual : Result, expected : Result, message : String) : Unit {
        EqualityFactR(actual, expected, message);
    }

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Diagnostics.AllEqualityFactB" instead.
    @Deprecated("Microsoft.Quantum.Diagnostics.AllEqualityFactB")
    function AssertBoolArrayEqual (actual : Bool[], expected : Bool[], message : String) : Unit {
        AllEqualityFactB(actual, expected, message);
    }

}
