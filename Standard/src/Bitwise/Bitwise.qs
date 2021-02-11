// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Bitwise {

    /// # Summary
    /// Shifts the bitwise representation of a number left by a given number of
    /// bits.
    ///
    /// # Input
    /// ## value
    /// The number whose bitwise representation is to be shifted to the left
    /// (more significant).
    /// ## amount
    /// The number of bits by which `value` is to be shifted to the left.
    ///
    /// # Output
    /// The value of `value`, shifted left by `amount` bits.
    ///
    /// # Remarks
    /// The following are equivalent:
    /// ```qsharp
    /// let c = a <<< b;
    /// let c = LeftShiftedI(a, b);
    /// ```
    function LeftShiftedI(value : Int, amount : Int) : Int {
        return value <<< amount;
    }

    /// # Summary
    /// Shifts the bitwise representation of a number left by a given number of
    /// bits.
    ///
    /// # Input
    /// ## value
    /// The number whose bitwise representation is to be shifted to the left
    /// (more significant).
    /// ## amount
    /// The number of bits by which `value` is to be shifted to the left.
    ///
    /// # Output
    /// The value of `value`, shifted left by `amount` bits.
    ///
    /// # Remarks
    /// The following are equivalent:
    /// ```qsharp
    /// let c = a <<< b;
    /// let c = LeftShiftedL(a, b);
    /// ```
    function LeftShiftedL(value : BigInt, amount : Int) : BigInt {
        return value <<< amount;
    }

    /// # Summary
    /// Shifts the bitwise representation of a number right by a given number of
    /// bits.
    ///
    /// # Input
    /// ## value
    /// The number whose bitwise representation is to be shifted to the right
    /// (less significant).
    /// ## amount
    /// The number of bits by which `value` is to be shifted to the right.
    ///
    /// # Output
    /// The value of `value`, shifted right by `amount` bits.
    ///
    /// # Remarks
    /// The following are equivalent:
    /// ```qsharp
    /// let c = a >>> b;
    /// let c = RightShiftedI(a, b);
    /// ```
    function RightShiftedI(value : Int, amount : Int) : Int {
        return value >>> amount;
    }

    /// # Summary
    /// Shifts the bitwise representation of a number right by a given number of
    /// bits.
    ///
    /// # Input
    /// ## value
    /// The number whose bitwise representation is to be shifted to the right
    /// (less significant).
    /// ## amount
    /// The number of bits by which `value` is to be shifted to the right.
    ///
    /// # Output
    /// The value of `value`, shifted right by `amount` bits.
    ///
    /// # Remarks
    /// The following are equivalent:
    /// ```qsharp
    /// let c = a >>> b;
    /// let c = RightShiftedL(a, b);
    /// ```
    function RightShiftedL(value : BigInt, amount : Int) : BigInt {
        return value >>> amount;
    }
}
