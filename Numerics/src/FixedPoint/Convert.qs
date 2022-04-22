// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Convert {
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Math;

    /// # Summary
    /// Computes fixed-point approximation for a double and returns it as `Bool` array.
    ///
    /// # Input
    /// ## integerBits
    /// Assumed number of integer bits (including the sign bit).
    /// ## fractionalBits
    /// Assumed number of fractional bits.
    /// ## value
    /// Value to be approximated.
    function FixedPointAsBoolArray(integerBits : Int, fractionalBits : Int, value : Double) : Bool[] {
        let numBits = integerBits + fractionalBits;
        let sign = value < 0.0;

        mutable result = [false, size = numBits];
        mutable rescaledConstant = PowD(2.0, IntAsDouble(fractionalBits)) * AbsD(value) + 0.5;
        mutable keepAdding = sign;

        for idx in 0..numBits - 1 {
            let intConstant = Floor(rescaledConstant);
            set rescaledConstant = rescaledConstant / 2.0;
            mutable currentBit = (intConstant &&& 1) == (sign ? 0 | 1);
            if keepAdding {
                set keepAdding = currentBit;
                set currentBit = not currentBit;
            }
            if currentBit {
                set result w/= idx <- true;
            }
        }

        return result;
    }

    /// # Summary
    /// Returns the double value of a fixed-point approximation from of a `Bool` array.
    ///
    /// # Input
    /// ## integerBits
    /// Assumed number of integer bits (including the sign bit).
    /// ## bits
    /// Bit-string representation of approximated number.
    function BoolArrayAsFixedPoint(integerBits : Int, bits : Bool[]) : Double {
        let numBits = Length(bits);
        let intPart = (Tail(bits) ? -(1 <<< (numBits - 1)) | 0) + BoolArrayAsInt(Most(bits));
        return IntAsDouble(intPart) / PowD(2.0, IntAsDouble(numBits - integerBits));
    }

    /// # Summary
    /// Discretizes a double value as a fixed-point approximation and returns its value as a double.
    ///
    /// # Input
    /// ## integerBits
    /// Assumed number of integer bits (including the sign bit).
    /// ## fractionalBits
    /// Assumed number of fractional bits.
    /// ## value
    /// Value to be approximated.
    function DoubleAsFixedPoint(integerBits : Int, fractionalBits : Int, value : Double) : Double {
        return BoolArrayAsFixedPoint(integerBits, FixedPointAsBoolArray(integerBits, fractionalBits, value));
    }

}
