// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Math {

    /// # Summary
    /// Returns the unary negation of an input.
    ///
    /// # Input
    /// ## input
    /// A value whose negation is to be returned.
    ///
    /// # Output
    /// The unary negation of `input`.
    function NegationI(input : Int) : Int {
        return -input;
    }

    /// # Summary
    /// Returns the unary negation of an input.
    ///
    /// # Input
    /// ## input
    /// A value whose negation is to be returned.
    ///
    /// # Output
    /// The unary negation of `input`.
    function NegationD(input : Double) : Double {
        return -input;
    }

    /// # Summary
    /// Returns the unary negation of an input.
    ///
    /// # Input
    /// ## input
    /// A value whose negation is to be returned.
    ///
    /// # Output
    /// The unary negation of `input`.
    function NegationL(input : BigInt) : BigInt {
        return -input;
    }

    /// # Summary
    /// Returns the sum of two inputs.
    ///
    /// # Input
    /// ## a
    /// The first input $a$ to be summed.
    /// ## b
    /// The second input $b$ to be summed.
    ///
    /// # Output
    /// The sum $a + b$.
    function PlusI(a : Int, b : Int) : Int {
        return a + b;
    }

    /// # Summary
    /// Returns the sum of two inputs.
    ///
    /// # Input
    /// ## a
    /// The first input $a$ to be summed.
    /// ## b
    /// The second input $b$ to be summed.
    ///
    /// # Output
    /// The sum $a + b$.
    function PlusD(a : Double, b : Double) : Double {
        return a + b;
    }

    /// # Summary
    /// Returns the sum of two inputs.
    ///
    /// # Input
    /// ## a
    /// The first input $a$ to be summed.
    /// ## b
    /// The second input $b$ to be summed.
    ///
    /// # Output
    /// The sum $a + b$.
    function PlusL(a : BigInt, b : BigInt) : BigInt {
        return a + b;
    }

    /// # Summary
    /// Returns the sum of two inputs.
    ///
    /// # Input
    /// ## a
    /// The first input $a$ to be summed.
    /// ## b
    /// The second input $b$ to be summed.
    ///
    /// # Output
    /// The sum $a + b$.
    function PlusA<'Element>(a : 'Element[], b : 'Element[]) : 'Element[] {
        return a + b;
    }

    /// # Summary
    /// Returns the difference between two inputs.
    ///
    /// # Input
    /// ## a
    /// The first input $a$ to be subtracted.
    /// ## b
    /// The second input $b$ to be subtracted.
    ///
    /// # Output
    /// The difference $a - b$.
    function MinusI(a : Int, b : Int) : Int {
        return a - b;
    }

    /// # Summary
    /// Returns the difference between two inputs.
    ///
    /// # Input
    /// ## a
    /// The first input $a$ to be subtracted.
    /// ## b
    /// The second input $b$ to be subtracted.
    ///
    /// # Output
    /// The difference $a - b$.
    function MinusD(a : Double, b : Double) : Double {
        return a - b;
    }

    /// # Summary
    /// Returns the difference between two inputs.
    ///
    /// # Input
    /// ## a
    /// The first input $a$ to be subtracted.
    /// ## b
    /// The second input $b$ to be subtracted.
    ///
    /// # Output
    /// The difference $a - b$.
    function MinusL(a : BigInt, b : BigInt) : BigInt {
        return a - b;
    }

    /// # Summary
    /// Returns the product of two inputs.
    ///
    /// # Input
    /// ## a
    /// The first input $a$ to be multiplied.
    /// ## b
    /// The second input $b$ to be multiplied.
    ///
    /// # Output
    /// The product $a \times b$.
    function TimesI(a : Int, b : Int) : Int {
        return a * b;
    }

    /// # Summary
    /// Returns the product of two inputs.
    ///
    /// # Input
    /// ## a
    /// The first input $a$ to be multiplied.
    /// ## b
    /// The second input $b$ to be multiplied.
    ///
    /// # Output
    /// The product $a \times b$.
    function TimesD(a : Double, b : Double) : Double {
        return a * b;
    }

    /// # Summary
    /// Returns the product of two inputs.
    ///
    /// # Input
    /// ## a
    /// The first input $a$ to be multiplied.
    /// ## b
    /// The second input $b$ to be multiplied.
    ///
    /// # Output
    /// The product $a \times b$.
    function TimesL(a : BigInt, b : BigInt) : BigInt {
        return a * b;
    }

    /// # Summary
    /// Returns the quotient of two inputs.
    ///
    /// # Input
    /// ## a
    /// The first input $a$ to be divided.
    /// ## b
    /// The second input $b$ to be divided.
    ///
    /// # Output
    /// The quotient $a / b$.
    function DividedByI(a : Int, b : Int) : Int {
        return a / b;
    }

    /// # Summary
    /// Returns the quotient of two inputs.
    ///
    /// # Input
    /// ## a
    /// The first input $a$ to be divided.
    /// ## b
    /// The second input $b$ to be divided.
    ///
    /// # Output
    /// The quotient $a / b$.
    function DividedByD(a : Double, b : Double) : Double {
        return a / b;
    }

    /// # Summary
    /// Returns the quotient of two inputs.
    ///
    /// # Input
    /// ## a
    /// The first input $a$ to be divided.
    /// ## b
    /// The second input $b$ to be divided.
    ///
    /// # Output
    /// The quotient $a / b$.
    function DividedByL(a : BigInt, b : BigInt) : BigInt {
        return a / b;
    }

    /// # Summary
    /// Returns the quotient of two inputs.
    ///
    /// # Input
    /// ## a
    /// The first input $a$ to be divided.
    /// ## b
    /// The second input $b$ to be divided.
    ///
    /// # Output
    /// The quotient $a / b$.
    function DividedByC(a : Complex, b : Complex) : Complex {
        return TimesC(a, PowC(b, -1.0));
    }

    /// # Summary
    /// Returns the quotient of two inputs.
    ///
    /// # Input
    /// ## a
    /// The first input $a$ to be divided.
    /// ## b
    /// The second input $b$ to be divided.
    ///
    /// # Output
    /// The quotient $a / b$.
    function DividedByCP(a : ComplexPolar, b : ComplexPolar) : ComplexPolar {
        return TimesCP(a, PowCP(b, -1.0));
    }

    /// # Summary
    /// Returns the modulus of a number with respect to another number.
    ///
    /// # Input
    /// ## a
    /// The input $a$ whose modulus is to be returned.
    /// ## b
    /// The number with respect to which the modulus of $a$ is to be returned.
    ///
    /// # Output
    /// The modulus $a \operatorname{mod} b$.
    function ModI(a : Int, b : Int) : Int {
        return a % b;
    }

    /// # Summary
    /// Returns the modulus of a number with respect to another number.
    ///
    /// # Input
    /// ## a
    /// The input $a$ whose modulus is to be returned.
    /// ## b
    /// The number with respect to which the modulus of $a$ is to be returned.
    ///
    /// # Output
    /// The modulus $a \operatorname{mod} b$.
    function ModL(a : BigInt, b : BigInt) : BigInt {
        return a % b;
    }

    /// # Summary
    /// Returns a number raised to a given power.
    ///
    /// # Input
    /// ## a
    /// The number $a$ that is to be raised.
    /// ## power
    /// The power $b$ to which $a$ should be raised.
    ///
    /// # Output
    /// The power $a^b$
    function PowI(a : Int, power : Int) : Int {
        return a ^ power;
    }

    /// # Summary
    /// Returns a number raised to a given power.
    ///
    /// # Input
    /// ## a
    /// The number $a$ that is to be raised.
    /// ## power
    /// The power $b$ to which $a$ should be raised.
    ///
    /// # Output
    /// The power $a^b$
    function PowL(a : BigInt, power : Int) : BigInt {
        return a ^ power;
    }

}
