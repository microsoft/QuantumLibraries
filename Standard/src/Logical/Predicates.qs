// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Logical {
    open Microsoft.Quantum.Math;

    /// # Summary
    /// Returns true if and only if two inputs are equal.
    ///
    /// # Input
    /// ## a
    /// The first value to be compared.
    /// ## b
    /// The second value to be compared.
    ///
    /// # Output
    /// `true` if and only if `a` is equal to `b`.
    ///
    /// # Remarks
    /// The following are equivalent:
    /// ```Q#
    /// let cond = a == b;
    /// let cond = EqualI(a, b);
    /// ```
    function EqualI(a : Int, b : Int) : Bool {
        return a == b;
    }

    /// # Summary
    /// Returns true if and only if two inputs are equal.
    ///
    /// # Input
    /// ## a
    /// The first value to be compared.
    /// ## b
    /// The second value to be compared.
    ///
    /// # Output
    /// `true` if and only if `a` is equal to `b`.
    ///
    /// # Remarks
    /// The following are equivalent:
    /// ```Q#
    /// let cond = a == b;
    /// let cond = EqualL(a, b);
    /// ```
    function EqualL(a : BigInt, b : BigInt) : Bool {
        return a == b;
    }

    /// # Summary
    /// Returns true if and only if two inputs are equal.
    ///
    /// # Input
    /// ## a
    /// The first value to be compared.
    /// ## b
    /// The second value to be compared.
    ///
    /// # Output
    /// `true` if and only if `a` is equal to `b`.
    ///
    /// # Remarks
    /// The following are equivalent:
    /// ```Q#
    /// let cond = a == b;
    /// let cond = EqualD(a, b);
    /// ```
    function EqualD(a : Double, b : Double) : Bool {
        return a == b;
    }

    /// # Summary
    /// Returns true if and only if two inputs are nearly equal (that is, within
    /// a tolerance of 1e-12).
    ///
    /// # Input
    /// ## a
    /// The first value to be compared.
    /// ## b
    /// The second value to be compared.
    ///
    /// # Output
    /// `true` if and only if `a` is nearly equal to `b`.
    ///
    /// # Remarks
    /// The following are equivalent:
    /// ```Q#
    /// let cond = Microsoft.Quantum.Math.AbsD(a - b) < 1e-12;
    /// let cond = NearlyEqualD(a, b);
    /// ```
    function NearlyEqualD(a : Double, b : Double) : Bool {
        return AbsD(a - b) < 1e-12;
    }

    /// # Summary
    /// Returns true if and only if two inputs are equal.
    ///
    /// # Input
    /// ## a
    /// The first value to be compared.
    /// ## b
    /// The second value to be compared.
    ///
    /// # Output
    /// `true` if and only if `a` is equal to `b`.
    ///
    /// # Remarks
    /// The following are equivalent:
    /// ```Q#
    /// let cond = a == b;
    /// let cond = EqualR(a, b);
    /// ```
    function EqualR(a : Result, b : Result) : Bool {
        return a == b;
    }

    /// # Summary
    /// Returns true if and only if two inputs are equal.
    ///
    /// # Input
    /// ## a
    /// The first value to be compared.
    /// ## b
    /// The second value to be compared.
    ///
    /// # Output
    /// `true` if and only if `a` is equal to `b`.
    ///
    /// # Remarks
    /// The following are equivalent:
    /// ```Q#
    /// let cond = a == b;
    /// let cond = EqualB(a, b);
    /// ```
    function EqualB(a : Bool, b : Bool) : Bool {
        return a == b;
    }

    /// # Summary
    /// Returns true if and only if two inputs are equal.
    ///
    /// # Input
    /// ## a
    /// The first value to be compared.
    /// ## b
    /// The second value to be compared.
    ///
    /// # Output
    /// `true` if and only if `a` is equal to `b`.
    function EqualC(a : Complex, b : Complex) : Bool {
        let ((reA, imA), (reB, imB)) = (a!, b!);
        return reA == reB and imA == imB;
    }

    /// # Summary
    /// Returns true if and only if two inputs are equal.
    ///
    /// # Input
    /// ## a
    /// The first value to be compared.
    /// ## b
    /// The second value to be compared.
    ///
    /// # Output
    /// `true` if and only if `a` is equal to `b`.
    function EqualCP(a : ComplexPolar, b : ComplexPolar) : Bool {
        return EqualC(ComplexPolarAsComplex(a), ComplexPolarAsComplex(b));
    }

    /// # Summary
    /// Returns true if and only if two inputs are not equal.
    ///
    /// # Input
    /// ## a
    /// The first value to be compared.
    /// ## b
    /// The second value to be compared.
    ///
    /// # Output
    /// `true` if and only if `a` is not equal to `b`.
    ///
    /// # Remarks
    /// The following are equivalent:
    /// ```Q#
    /// let cond = a != b;
    /// let cond = NotEqualI(a, b);
    /// ```
    function NotEqualI(a : Int, b : Int) : Bool {
        return a != b;
    }

    /// # Summary
    /// Returns true if and only if two inputs are not equal.
    ///
    /// # Input
    /// ## a
    /// The first value to be compared.
    /// ## b
    /// The second value to be compared.
    ///
    /// # Output
    /// `true` if and only if `a` is not equal to `b`.
    ///
    /// # Remarks
    /// The following are equivalent:
    /// ```Q#
    /// let cond = a != b;
    /// let cond = NotEqualL(a, b);
    /// ```
    function NotEqualL(a : BigInt, b : BigInt) : Bool {
        return a != b;
    }

    /// # Summary
    /// Returns true if and only if two inputs are not equal.
    ///
    /// # Input
    /// ## a
    /// The first value to be compared.
    /// ## b
    /// The second value to be compared.
    ///
    /// # Output
    /// `true` if and only if `a` is not equal to `b`.
    ///
    /// # Remarks
    /// The following are equivalent:
    /// ```Q#
    /// let cond = a != b;
    /// let cond = NotEqualD(a, b);
    /// ```
    function NotEqualD(a : Double, b : Double) : Bool {
        return a != b;
    }

    /// # Summary
    /// Returns true if and only if two inputs are not nearly equal (that is,
    /// are not within a tolerance of 1e-12).
    ///
    /// # Input
    /// ## a
    /// The first value to be compared.
    /// ## b
    /// The second value to be compared.
    ///
    /// # Output
    /// `true` if and only if `a` is not nearly equal to `b`.
    ///
    /// # Remarks
    /// The following are equivalent:
    /// ```Q#
    /// let cond = Microsoft.Quantum.Math.AbsD(a - b) >= 1e-12;
    /// let cond = NotNearlyEqualD(a, b);
    /// ```
    function NotNearlyEqualD(a : Double, b : Double) : Bool {
        return not NearlyEqualD(a, b);
    }

    /// # Summary
    /// Returns true if and only if two inputs are not equal.
    ///
    /// # Input
    /// ## a
    /// The first value to be compared.
    /// ## b
    /// The second value to be compared.
    ///
    /// # Output
    /// `true` if and only if `a` is not equal to `b`.
    ///
    /// # Remarks
    /// The following are equivalent:
    /// ```Q#
    /// let cond = a != b;
    /// let cond = NotEqualR(a, b);
    /// ```
    function NotEqualR(a : Result, b : Result) : Bool {
        return a != b;
    }

    /// # Summary
    /// Returns true if and only if two inputs are not equal.
    ///
    /// # Input
    /// ## a
    /// The first value to be compared.
    /// ## b
    /// The second value to be compared.
    ///
    /// # Output
    /// `true` if and only if `a` is not equal to `b`.
    ///
    /// # Remarks
    /// The following are equivalent:
    /// ```Q#
    /// let cond = a != b;
    /// let cond = NotEqualB(a, b);
    /// ```
    function NotEqualB(a : Bool, b : Bool) : Bool {
        return a != b;
    }

    /// # Summary
    /// Returns true if and only if two inputs are not equal.
    ///
    /// # Input
    /// ## a
    /// The first value to be compared.
    /// ## b
    /// The second value to be compared.
    ///
    /// # Output
    /// `true` if and only if `a` is not equal to `b`.
    function NotEqualC(a : Complex, b : Complex) : Bool {
        return not EqualC(a, b);
    }

    /// # Summary
    /// Returns true if and only if two inputs are not equal.
    ///
    /// # Input
    /// ## a
    /// The first value to be compared.
    /// ## b
    /// The second value to be compared.
    ///
    /// # Output
    /// `true` if and only if `a` is not equal to `b`.
    function NotEqualCP(a : ComplexPolar, b : ComplexPolar) : Bool {
        return not EqualCP(a, b);
    }

    /// # Summary
    /// Returns true if and only if a number is greater than another number.
    ///
    /// # Input
    /// ## a
    /// The first value to be compared.
    /// ## b
    /// The second value to be compared.
    ///
    /// # Output
    /// `true` if and only if `a` is strictly greater than `b`.
    ///
    /// # Remarks
    /// The following are equivalent:
    /// ```Q#
    /// let cond = a > b;
    /// let cond = GreaterThanI(a, b);
    /// ```
    function GreaterThanI(a : Int, b : Int) : Bool {
        return a > b;
    }

    /// # Summary
    /// Returns true if and only if a number is greater than another number.
    ///
    /// # Input
    /// ## a
    /// The first value to be compared.
    /// ## b
    /// The second value to be compared.
    ///
    /// # Output
    /// `true` if and only if `a` is strictly greater than `b`.
    ///
    /// # Remarks
    /// The following are equivalent:
    /// ```Q#
    /// let cond = a > b;
    /// let cond = GreaterThanL(a, b);
    /// ```
    function GreaterThanL(a : BigInt, b : BigInt) : Bool {
        return a > b;
    }

    /// # Summary
    /// Returns true if and only if a number is greater than another number.
    ///
    /// # Input
    /// ## a
    /// The first value to be compared.
    /// ## b
    /// The second value to be compared.
    ///
    /// # Output
    /// `true` if and only if `a` is strictly greater than `b`.
    ///
    /// # Remarks
    /// The following are equivalent:
    /// ```Q#
    /// let cond = a > b;
    /// let cond = GreaterThanD(a, b);
    /// ```
    function GreaterThanD(a : Double, b : Double) : Bool {
        return a > b;
    }

    /// # Summary
    /// Returns true if and only if a number is greater than or equal to another
    /// number.
    ///
    /// # Input
    /// ## a
    /// The first value to be compared.
    /// ## b
    /// The second value to be compared.
    ///
    /// # Output
    /// `true` if and only if `a` is greater than or is equal to `b`.
    ///
    /// # Remarks
    /// The following are equivalent:
    /// ```Q#
    /// let cond = a >= b;
    /// let cond = GreaterThanOrEqualI(a, b);
    /// ```
    function GreaterThanOrEqualI(a : Int, b : Int) : Bool {
        return a >= b;
    }

    /// # Summary
    /// Returns true if and only if a number is greater than or equal to another
    /// number.
    ///
    /// # Input
    /// ## a
    /// The first value to be compared.
    /// ## b
    /// The second value to be compared.
    ///
    /// # Output
    /// `true` if and only if `a` is greater than or is equal to `b`.
    ///
    /// # Remarks
    /// The following are equivalent:
    /// ```Q#
    /// let cond = a >= b;
    /// let cond = GreaterThanOrEqualL(a, b);
    /// ```
    function GreaterThanOrEqualL(a : BigInt, b : BigInt) : Bool {
        return a >= b;
    }

    /// # Summary
    /// Returns true if and only if a number is greater than or equal to another
    /// number.
    ///
    /// # Input
    /// ## a
    /// The first value to be compared.
    /// ## b
    /// The second value to be compared.
    ///
    /// # Output
    /// `true` if and only if `a` is greater than or is equal to `b`.
    ///
    /// # Remarks
    /// The following are equivalent:
    /// ```Q#
    /// let cond = a >= b;
    /// let cond = GreaterThanOrEqualD(a, b);
    /// ```
    function GreaterThanOrEqualD(a : Double, b : Double) : Bool {
        return a >= b;
    }

    /// # Summary
    /// Returns true if and only if a number is less than another number.
    ///
    /// # Input
    /// ## a
    /// The first value to be compared.
    /// ## b
    /// The second value to be compared.
    ///
    /// # Output
    /// `true` if and only if `a` is strictly less than `b`.
    ///
    /// # Remarks
    /// The following are equivalent:
    /// ```Q#
    /// let cond = a < b;
    /// let cond = LessThanI(a, b);
    /// ```
    function LessThanI(a : Int, b : Int) : Bool {
        return a < b;
    }

    /// # Summary
    /// Returns true if and only if a number is less than another number.
    ///
    /// # Input
    /// ## a
    /// The first value to be compared.
    /// ## b
    /// The second value to be compared.
    ///
    /// # Output
    /// `true` if and only if `a` is strictly less than `b`.
    ///
    /// # Remarks
    /// The following are equivalent:
    /// ```Q#
    /// let cond = a < b;
    /// let cond = LessThanL(a, b);
    /// ```
    function LessThanL(a : BigInt, b : BigInt) : Bool {
        return a < b;
    }

    /// # Summary
    /// Returns true if and only if a number is less than another number.
    ///
    /// # Input
    /// ## a
    /// The first value to be compared.
    /// ## b
    /// The second value to be compared.
    ///
    /// # Output
    /// `true` if and only if `a` is strictly less than `b`.
    ///
    /// # Remarks
    /// The following are equivalent:
    /// ```Q#
    /// let cond = a < b;
    /// let cond = LessThanD(a, b);
    /// ```
    function LessThanD(a : Double, b : Double) : Bool {
        return a < b;
    }

    /// # Summary
    /// Returns true if and only if a number is less than or equal to another
    /// number.
    ///
    /// # Input
    /// ## a
    /// The first value to be compared.
    /// ## b
    /// The second value to be compared.
    ///
    /// # Output
    /// `true` if and only if `a` is less than or equal to `b`.
    ///
    /// # Remarks
    /// The following are equivalent:
    /// ```Q#
    /// let cond = a <= b;
    /// let cond = LessThanOrEqualI(a, b);
    /// ```
    function LessThanOrEqualI(a : Int, b : Int) : Bool {
        return a <= b;
    }

    /// # Summary
    /// Returns true if and only if a number is less than or equal to another
    /// number.
    ///
    /// # Input
    /// ## a
    /// The first value to be compared.
    /// ## b
    /// The second value to be compared.
    ///
    /// # Output
    /// `true` if and only if `a` is less than or equal to `b`.
    ///
    /// # Remarks
    /// The following are equivalent:
    /// ```Q#
    /// let cond = a <= b;
    /// let cond = LessThanOrEqualL(a, b);
    /// ```
    function LessThanOrEqualL(a : BigInt, b : BigInt) : Bool {
        return a <= b;
    }

    /// # Summary
    /// Returns true if and only if a number is less than or equal to another
    /// number.
    ///
    /// # Input
    /// ## a
    /// The first value to be compared.
    /// ## b
    /// The second value to be compared.
    ///
    /// # Output
    /// `true` if and only if `a` is less than or equal to `b`.
    ///
    /// # Remarks
    /// The following are equivalent:
    /// ```Q#
    /// let cond = a <= b;
    /// let cond = LessThanOrEqualD(a, b);
    /// ```
    function LessThanOrEqualD(a : Double, b : Double) : Bool {
        return a <= b;
    }

}
