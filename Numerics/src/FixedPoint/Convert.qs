// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Convert {
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Math;

    /// # Summary
    /// Returns the double value of a fixed-point approximation from of a `Bool` array.
    ///
    /// # Input
    /// ## integerBits
    /// Assumed number of integerBits (including the sign big)
    /// ## bits
    /// Bit-string representation of approximated number
    internal function BoolArrayAsFixedPoint(integerBits : Int, bits : Bool[]) : Double {
        let numBits = Length(bits);
        let intPart = (Tail(bits) ? -(1 <<< (numBits - 1)) | 0) + BoolArrayAsInt(Most(bits));
        return IntAsDouble(intPart) / PowD(2.0, IntAsDouble(numBits - integerBits));
    }

}
