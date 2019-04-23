// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {
    open Microsoft.Quantum.Math;
    open Microsoft.Quantum.Warnings;

    /// # Deprecated
    /// This function has been removed.
    function IntAbs(input : Int) : Int {
        // TODO: add call to _Removed.
        return AbsI(input);
    }

    /// # Deprecated
    /// This function has been removed.
    function IntMax(a : Int, b : Int) : Int {
        // TODO: add call to _Removed.
        return MaxI(a, b);
    }

    /// # Deprecated
    /// Please use @"microsoft.quantum.math.pnormalized".
    function PNormalize(p : Double, array : Double[]) : Double[] {
        _Renamed("Microsoft.Quantum.Canon.PNormalize", "Microsoft.Quantum.Math.PNormalized");
        return PNormalized(p, array);
    }
}
