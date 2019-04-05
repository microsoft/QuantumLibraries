// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {
    open Microsoft.Quantum.Primitive;
    open Microsoft.Quantum.Diagnostics;

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Diagnostics.EqualityWithinToleranceFact" instead.
    function AssertAlmostEqualTol (actual : Double, expected : Double, tolerance : Double) : Unit {
        Renamed("Microsoft.Quantum.Canon.AssertAlmostEqualTol", "Microsoft.Quantum.Diagnostics.EqualityWithinToleranceFact");
        return EqualityWithinToleranceFact(actual, expected, tolerance);
    }

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Diagnostics.NearEqualityFact" instead.
    function AssertAlmostEqual(actual : Double, expected : Double) : Unit {
        Renamed("Microsoft.Quantum.Canon.AssertAlmostEqual", "Microsoft.Quantum.Diagnostics.NearEqualityFact");
        NearEqualityFact(actual, expected);
    }

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Diagnostics.EqualityFactI" instead.
    function AssertIntEqual (actual : Int, expected : Int, message : String) : Unit {
        Renamed("Microsoft.Quantum.Canon.AssertIntEqual", "Microsoft.Quantum.Diagnostics.EqualityFactI");
        EqualityFactI(actual, expected, message);
    }

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Diagnostics.EqualityFactB" instead.
    function AssertBoolEqual(actual : Bool, expected : Bool, message : String) : Unit {
        Renamed("Microsoft.Quantum.Canon.AssertBoolEqual", "Microsoft.Quantum.Diagnostics.EqualityFactB");
        EqualityFactB(actual, expected, message);
    }

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Diagnostics.EqualityFactR" instead.
    function AssertResultEqual(actual : Result, expected : Result, message : String) : Unit {
        Renamed("Microsoft.Quantum.Canon.AssertResultEqual", "Microsoft.Quantum.Diagnostics.EqualityFactR");
        EqualityFactR(actual, expected, message);
    }

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Diagnostics.AllEqualityFactB" instead.
    function AssertBoolArrayEqual (actual : Bool[], expected : Bool[], message : String) : Unit {
        Renamed("Microsoft.Quantum.Canon.AssertBoolArrayEqual", "Microsoft.Quantum.Diagnostics.AllEqualityFactB");
        AllEqualityFactB(actual, expected, message);
    }

}
