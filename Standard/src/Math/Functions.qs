// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Math {
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Convert;
    open Microsoft.Quantum.Diagnostics;
    open Microsoft.Quantum.Arrays;

    /// # Summary
    /// Computes the base-2 logarithm of a number.
    ///
    /// # Input
    /// ## input
    /// A real number $x$.
    ///
    /// # Output
    /// The base-2 logarithm $y = \log_2(x)$ such that $x = 2^y$.
    function Lg (input : Double) : Double {
        return Log(input) / LogOf2();
    }

    /// # Summary
    /// Given an array of integers, returns the largest element.
    ///
    /// # Input
    /// ## values
    /// An array to take the maximum of.
    ///
    /// # Output
    /// The largest element of `values`.
    function Max (values : Int[]) : Int {
        mutable max = values[0];
        let nTerms = Length(values);

        for idx in 0 .. nTerms - 1 {
            if values[idx] > max {
                set max = values[idx];
            }
        }

        return max;
    }

    /// # Summary
    /// Given an array of integers, returns the smallest element.
    ///
    /// # Input
    /// ## values
    /// An array to take the minimum of.
    ///
    /// # Output
    /// The smallest element of `values`.
    function Min (values : Int[]) : Int {
        mutable min = values[0];
        let nTerms = Length(values);

        for idx in 0 .. nTerms - 1 {
            if values[idx] < min {
                set min = values[idx];
            }
        }

        return min;
    }


    /// # Summary
    /// Computes the modulus between two real numbers.
    ///
    /// # Input
    /// ## value
    /// A real number $x$ to take the modulus of.
    /// ## modulo
    /// A real number to take the modulus of $x$ with respect to.
    /// ## minValue
    /// The smallest value to be returned by this function.
    ///
    /// # Remarks
    /// This function computes the real modulus by wrapping the real
    /// line about the unit circle, then finding the angle on the
    /// unit circle corresponding to the input.
    /// The `minValue` input then effectively specifies where to cut the
    /// unit circle.
    ///
    /// ## Example
    /// ```qsharp
    ///     // Returns 3 π / 2.
    ///     let y = RealMod(5.5 * PI(), 2.0 * PI(), 0.0);
    ///     // Returns -1.2, since +3.6 and -1.2 are 4.8 apart on the real line,
    ///     // which is a multiple of 2.4.
    ///     let z = RealMod(3.6, 2.4, -1.2);
    /// ```
    function RealMod(value : Double, modulo : Double, minValue : Double) : Double
    {
        let fractionalValue = (2.0 * PI()) * ((value - minValue) / modulo - 0.5);
        let cosFracValue = Cos(fractionalValue);
        let sinFracValue = Sin(fractionalValue);
        let moduloValue = 0.5 + ArcTan2(sinFracValue, cosFracValue) / (2.0 * PI());
        let output = moduloValue * modulo + minValue;
        return output;
    }


    // NB: .NET's Math library does not provide hyperbolic arcfunctions.

    /// # Summary
    /// Computes the inverse hyperbolic cosine of a number.
    ///
    /// # Input
    /// ## x
    /// A real number $x\geq 1$.
    ///
    /// # Output
    /// A real number $y$ such that $x = \cosh(y)$.
    function ArcCosh (x : Double) : Double {
        return Log(x + Sqrt(x * x - 1.0));
    }


    /// # Summary
    /// Computes the inverse hyperbolic sine of a number.
    ///
    /// # Input
    /// ## x
    /// A real number $x$.
    ///
    /// # Output
    /// A real number $y$ such that $x = \operatorname{sinh}(y)$.
    function ArcSinh (x : Double) : Double
    {
        return Log(x + Sqrt(x * x + 1.0));
    }


    /// # Summary
    /// Computes the inverse hyperbolic tangent of a number.
    ///
    /// # Input
    /// ## x
    /// A real number $x$.
    ///
    /// # Output
    /// A real number $y$ such that $x = \tanh(y)$.
    function ArcTanh (x : Double) : Double
    {
        return Log((1.0 + x) / (1.0 - x)) * 0.5;
    }


    /// # Summary
    /// Computes the canonical residue of `value` modulo `modulus`.
    /// # Input
    /// ## value
    /// The value of which residue is computed
    /// ## modulus
    /// The modulus by which residues are take, must be positive
    /// # Output
    /// Integer $r$ between 0 and `modulus - 1` such that `value - r` is divisible by modulus
    ///
    /// # Remarks
    /// This function behaves different to how the operator `%` behaves in C# and Q# as in the result
    /// is always a non-negative integer between 0 and `modulus - 1`, even if value is negative.
    function ModulusI(value : Int, modulus : Int) : Int {
        Fact(modulus > 0, $"`modulus` must be positive");
        let r = value % modulus;
        return (r < 0) ? (r + modulus) | r;
    }

    /// # Summary
    /// Computes the canonical residue of `value` modulo `modulus`.
    /// # Input
    /// ## value
    /// The value of which residue is computed
    /// ## modulus
    /// The modulus by which residues are take, must be positive
    /// # Output
    /// Integer $r$ between 0 and `modulus - 1` such that `value - r` is divisible by modulus
    ///
    /// # Remarks
    /// This function behaves different to how the operator `%` behaves in C# and Q# as in the result
    /// is always a non-negative integer between 0 and `modulus - 1`, even if value is negative.
    function ModulusL(value : BigInt, modulus : BigInt) : BigInt {
        Fact(modulus > 0L, $"`modulus` must be positive");
        let r = value % modulus;
        return (r < 0L) ? (r + modulus) | r;
    }


    /// # Summary
    /// Returns an integer raised to a given power, with respect to a given
    /// modulus.
    ///
    /// # Description
    /// Let us denote expBase by $x$, power by $p$ and modulus by $N$.
    /// The function returns $x^p \operatorname{mod} N$.
    ///
    /// We assume that $N$, $x$ are positive and power is non-negative.
    ///
    /// # Remarks
    /// Takes time proportional to the number of bits in `power`, not the `power` itself.
    function ExpModI(expBase : Int, power : Int, modulus : Int) : Int {
        Fact(power >= 0, $"`power` must be non-negative");
        Fact(modulus > 0, $"`modulus` must be positive");
        Fact(expBase > 0, $"`expBase` must be positive");
        mutable res = 1;
        mutable expPow2mod = expBase;

        // express p as bit-string pₙ … p₀
        let powerBitExpansion = IntAsBoolArray(power, BitSizeI(power));
        let expBaseMod = expBase % modulus;

        for bit in powerBitExpansion {
            if bit {
                // if bit pₖ is 1, multiply res by expBase^(2ᵏ) (mod `modulus`)
                set res = (res * expPow2mod) % modulus;
            }

            // update value of expBase^(2ᵏ) (mod `modulus`)
            set expPow2mod = (expPow2mod * expPow2mod) % modulus;
        }

        return res;
    }

    /// # Summary
    /// Returns an integer raised to a given power, with respect to a given
    /// modulus.
    ///
    /// # Description
    /// Let us denote expBase by $x$, power by $p$ and modulus by $N$.
    /// The function returns $x^p \operatorname{mod} N$.
    ///
    /// We assume that $N$, $x$ are positive and power is non-negative.
    ///
    /// # Remarks
    /// Takes time proportional to the number of bits in `power`, not the `power` itself.
    function ExpModL(expBase : BigInt, power : BigInt, modulus : BigInt) : BigInt {
        Fact(power >= 0L, $"`power` must be non-negative");
        Fact(modulus > 0L, $"`modulus` must be positive");
        Fact(expBase > 0L, $"`expBase` must be positive");
        mutable res = 1L;
        mutable expPow2mod = expBase;

        // express p as bit-string pₙ … p₀
        let powerBitExpansion = BigIntAsBoolArray(power);
        let expBaseMod = expBase % modulus;

        for bit in powerBitExpansion {
            if bit {
                // if bit pₖ is 1, multiply res by expBase^(2ᵏ) (mod `modulus`)
                set res = (res * expPow2mod) % modulus;
            }

            // update value of expBase^(2ᵏ) (mod `modulus`)
            set expPow2mod = (expPow2mod * expPow2mod) % modulus;
        }

        return res;
    }

    /// # Summary
    /// Internal recursive call to calculate the GCD.
    function _ExtendedGreatestCommonDivisorI(signA : Int, signB : Int, r : (Int, Int), s : (Int, Int), t : (Int, Int)) : (Int, Int) {
        if Snd(r) == 0 {
            return (Fst(s) * signA, Fst(t) * signB);
        }

        let quotient = Fst(r) / Snd(r);
        let r_ = (Snd(r), Fst(r) - quotient * Snd(r));
        let s_ = (Snd(s), Fst(s) - quotient * Snd(s));
        let t_ = (Snd(t), Fst(t) - quotient * Snd(t));
        return _ExtendedGreatestCommonDivisorI(signA, signB, r_, s_, t_);
    }


    /// # Summary
    /// Returns the GCD of two integers, decomposed into a linear combination.
    ///
    /// # Description
    /// Computes a tuple $(u,v)$ such that $u \cdot a + v \cdot b = \operatorname{GCD}(a, b)$,
    /// where $\operatorname{GCD}$ is $a$
    /// greatest common divisor of $a$ and $b$. The GCD is always positive.
    ///
    /// # Input
    /// ## a
    /// the first number of which extended greatest common divisor is being computed
    /// ## b
    /// the second number of which extended greatest common divisor is being computed
    ///
    /// # Output
    /// Tuple $(u,v)$ with the property $u \cdot a + v \cdot b = \operatorname{GCD}(a, b)$.
    ///
    /// # References
    /// - This implementation is according to https://en.wikipedia.org/wiki/Extended_Euclidean_algorithm
    function ExtendedGreatestCommonDivisorI(a : Int, b : Int) : (Int, Int) {
        let signA = SignI(a);
        let signB = SignI(b);
        let s = (1, 0);
        let t = (0, 1);
        let r = (a * signA, b * signB);
        return _ExtendedGreatestCommonDivisorI(signA, signB, r, s, t);
    }

    /// # Summary
    /// Internal recursive call to calculate the GCD.
    function _ExtendedGreatestCommonDivisorL(signA : Int, signB : Int, r : (BigInt, BigInt), s : (BigInt, BigInt), t : (BigInt, BigInt)) : (BigInt, BigInt) {
        if Snd(r) == 0L {
            return (Fst(s) * IntAsBigInt(signA), Fst(t) * IntAsBigInt(signB));
        }

        let quotient = Fst(r) / Snd(r);
        let r_ = (Snd(r), Fst(r) - quotient * Snd(r));
        let s_ = (Snd(s), Fst(s) - quotient * Snd(s));
        let t_ = (Snd(t), Fst(t) - quotient * Snd(t));
        return _ExtendedGreatestCommonDivisorL(signA, signB, r_, s_, t_);
    }


    /// # Summary
    /// Returns the GCD of two integers, decomposed into a linear combination.
    ///
    /// # Description
    /// Computes a tuple $(u,v)$ such that $u \cdot a + v \cdot b = \operatorname{GCD}(a, b)$,
    /// where $\operatorname{GCD}$ is $a$
    /// greatest common divisor of $a$ and $b$. The GCD is always positive.
    ///
    /// # Input
    /// ## a
    /// the first number of which extended greatest common divisor is being computed
    /// ## b
    /// the second number of which extended greatest common divisor is being computed
    ///
    /// # Output
    /// Tuple $(u,v)$ with the property $u \cdot a + v \cdot b = \operatorname{GCD}(a, b)$.
    ///
    /// # References
    /// - This implementation is according to https://en.wikipedia.org/wiki/Extended_Euclidean_algorithm
    function ExtendedGreatestCommonDivisorL(a : BigInt, b : BigInt) : (BigInt, BigInt) {
        let signA = SignL(a);
        let signB = SignL(b);
        let s = (1l, 0L);
        let t = (0l, 1L);
        let r = (a * IntAsBigInt(signA), b * IntAsBigInt(signB));
        return _ExtendedGreatestCommonDivisorL(signA, signB, r, s, t);
    }


    /// # Summary
    /// Computes the greatest common divisor of two integers.
    ///
    /// # Description
    /// Computes the greatest common divisor of two integers $a$ and $b$.
    /// The GCD is always positive.
    ///
    /// # Input
    /// ## a
    /// the first number of which extended greatest common divisor is being computed
    /// ## b
    /// the second number of which extended greatest common divisor is being computed
    ///
    /// # Output
    /// Greatest common divisor of $a$ and $b$
    function GreatestCommonDivisorI(a : Int, b : Int) : Int {
        let (u, v) = ExtendedGreatestCommonDivisorI(a, b);
        return u * a + v * b;
    }

    /// # Summary
    /// Computes the greatest common divisor of two integers.
    ///
    /// # Description
    /// Computes the greatest common divisor of two integers $a$ and $b$.
    /// The GCD is always positive.
    ///
    /// # Input
    /// ## a
    /// the first number of which extended greatest common divisor is being computed
    /// ## b
    /// the second number of which extended greatest common divisor is being computed
    ///
    /// # Output
    /// Greatest common divisor of $a$ and $b$
    function GreatestCommonDivisorL(a : BigInt, b : BigInt) : BigInt {
        let (u, v) = ExtendedGreatestCommonDivisorL(a, b);
        return u * a + v * b;
    }


    /// # Summary
    /// Internal recursive call to calculate the GCD with a bound
    function _ContinuedFractionConvergentI(signA : Int, signB : Int, r : (Int, Int), s : (Int, Int), t : (Int, Int), denominatorBound : Int)
    : Fraction {
        if Snd(r) == 0 or AbsI(Snd(s)) > denominatorBound {
            return (Snd(r) == 0 and AbsI(Snd(s)) <= denominatorBound)
                   ? Fraction(-Snd(t) * signB, Snd(s) * signA)
                   | Fraction(-Fst(t) * signB, Fst(s) * signA);
        }

        let quotient = Fst(r) / Snd(r);
        let r_ = (Snd(r), Fst(r) - quotient * Snd(r));
        let s_ = (Snd(s), Fst(s) - quotient * Snd(s));
        let t_ = (Snd(t), Fst(t) - quotient * Snd(t));
        return _ContinuedFractionConvergentI(signA, signB, r_, s_, t_, denominatorBound);
    }


    /// # Summary
    /// Finds the continued fraction convergent closest to `fraction`
    /// with the denominator less or equal to `denominatorBound`
    ///
    /// # Input
    ///
    ///
    /// # Output
    /// Continued fraction closest to `fraction`
    /// with the denominator less or equal to `denominatorBound`
    function ContinuedFractionConvergentI(fraction : Fraction, denominatorBound : Int)
    : Fraction {
        Fact(denominatorBound > 0, $"Denominator bound must be positive");
        let (a, b) = fraction!;
        let signA = SignI(a);
        let signB = SignI(b);
        let s = (1, 0);
        let t = (0, 1);
        let r = (a * signA, b * signB);
        return _ContinuedFractionConvergentI(signA, signB, r, s, t, denominatorBound);
    }

    /// # Summary
    /// Internal recursive call to calculate the GCD with a bound
    function _ContinuedFractionConvergentL(signA : Int, signB : Int, r : (BigInt, BigInt), s : (BigInt, BigInt), t : (BigInt, BigInt), denominatorBound : BigInt) : BigFraction
    {
        if Snd(r) == 0L or AbsL(Snd(s)) > denominatorBound {
            return (Snd(r) == 0L and AbsL(Snd(s)) <= denominatorBound)
                   ? BigFraction(-Snd(t) * IntAsBigInt(signB), Snd(s) * IntAsBigInt(signA))
                   | BigFraction(-Fst(t) * IntAsBigInt(signB), Fst(s) * IntAsBigInt(signA));
        }

        let quotient = Fst(r) / Snd(r);
        let r_ = (Snd(r), Fst(r) - quotient * Snd(r));
        let s_ = (Snd(s), Fst(s) - quotient * Snd(s));
        let t_ = (Snd(t), Fst(t) - quotient * Snd(t));
        return _ContinuedFractionConvergentL(signA, signB, r_, s_, t_, denominatorBound);
    }


    /// # Summary
    /// Finds the continued fraction convergent closest to `fraction`
    /// with the denominator less or equal to `denominatorBound`
    ///
    /// # Input
    ///
    ///
    /// # Output
    /// Continued fraction closest to `fraction`
    /// with the denominator less or equal to `denominatorBound`
    function ContinuedFractionConvergentL(fraction : BigFraction, denominatorBound : BigInt)
    : BigFraction {
        Fact(denominatorBound > 0L, $"Denominator bound must be positive");
        let (a, b) = fraction!;
        let signA = SignL(a);
        let signB = SignL(b);
        let s = (1L, 0L);
        let t = (0L, 1L);
        let r = (a * IntAsBigInt(signA), b * IntAsBigInt(signB));
        return _ContinuedFractionConvergentL(signA, signB, r, s, t, denominatorBound);
    }

    /// # Summary
    /// Returns if two integers are co-prime.
    ///
    /// # Description
    /// Returns true if $a$ and $b$ are co-prime and false otherwise.
    ///
    /// # Input
    /// ## a
    /// the first number of which co-primality is being tested
    /// ## b
    /// the second number of which co-primality is being tested
    ///
    /// # Output
    /// True, if $a$ and $b$ are co-prime (e.g. their greatest common divisor is 1 ),
    /// and false otherwise
    function IsCoprimeI(a : Int, b : Int) : Bool {
        let (u, v) = ExtendedGreatestCommonDivisorI(a, b);
        return u * a + v * b == 1;
    }

    /// # Summary
    /// Returns if two integers are co-prime.
    ///
    /// # Description
    /// Returns true if $a$ and $b$ are co-prime and false otherwise.
    ///
    /// # Input
    /// ## a
    /// the first number of which co-primality is being tested
    /// ## b
    /// the second number of which co-primality is being tested
    ///
    /// # Output
    /// True, if $a$ and $b$ are co-prime (e.g. their greatest common divisor is 1 ),
    /// and false otherwise
    function IsCoprimeL(a : BigInt, b : BigInt) : Bool {
        let (u, v) = ExtendedGreatestCommonDivisorL(a, b);
        return u * a + v * b == 1L;
    }

    /// # Summary
    /// Returns the multiplicative inverse of a modular integer.
    ///
    /// # Description
    /// Returns $b$ such that $a \cdot b = 1 (\operatorname{mod} \texttt{modulus})$.
    ///
    /// # Input
    /// ## a
    /// The number being inverted
    /// ## modulus
    /// The modulus according to which the numbers are inverted
    ///
    /// # Output
    /// Integer $b$ such that $a \cdot b = 1 (\operatorname{mod} \texttt{modulus})$.
    function InverseModI(a : Int, modulus : Int) : Int {
        let (u, v) = ExtendedGreatestCommonDivisorI(a, modulus);
        let gcd = u * a + v * modulus;
        EqualityFactI(gcd, 1, $"`a` and `modulus` must be co-prime");
        return ModulusI(u, modulus);
    }

    /// # Summary
    /// Returns $b$ such that $a \cdot b = 1 (\operatorname{mod} \texttt{modulus})$.
    ///
    /// # Input
    /// ## a
    /// The number being inverted
    /// ## modulus
    /// The modulus according to which the numbers are inverted
    ///
    /// # Output
    /// Integer $b$ such that $a \cdot b = 1 (\operatorname{mod} \texttt{modulus})$.
    function InverseModL(a : BigInt, modulus : BigInt) : BigInt {
        let (u, v) = ExtendedGreatestCommonDivisorL(a, modulus);
        let gcd = u * a + v * modulus;
        EqualityFactL(gcd, 1L, $"`a` and `modulus` must be co-prime");
        return ModulusL(u, modulus);
    }


    /// # Summary
    /// Helper function used to recursively calculate the bitsize of a value.
    internal function AccumulatedBitsizeI(val : Int, bitsize : Int) : Int {
        return val == 0 ? bitsize | AccumulatedBitsizeI(val / 2, bitsize + 1);
    }


    /// # Summary
    /// For a non-negative integer `a`, returns the number of bits required to represent `a`.
    ///
    /// # Remarks
    /// This function returns the smallest $n$ such that $a < 2^n$.
    ///
    /// # Input
    /// ## a
    /// The integer whose bit-size is to be computed.
    ///
    /// # Output
    /// The bit-size of `a`.
    function BitSizeI(a : Int) : Int {
        Fact(a >= 0, $"`a` must be non-negative");
        return AccumulatedBitsizeI(a, 0);
    }


    /// # Summary
    /// For a non-negative integer `a`, returns the number of bits required to represent `a`.
    ///
    /// # Remarks
    /// This function returns the smallest $n$ such that $a < 2^n$.
    ///
    /// # Input
    /// ## a
    /// The integer whose bit-size is to be computed.
    ///
    /// # Output
    /// The bit-size of `a`.
    function BitSizeL(a : BigInt) : Int {
        Fact(a >= 0L, $"`a` must be non-negative");
        mutable bitsize = 0;
        mutable val = a;
        while val != 0L {
            set bitsize += 1;
            set val /= 2L;
        }
        return bitsize;
    }


    /// # Summary
    /// Returns the p-norm of a vector of real numbers.
    ///
    /// # Description
    /// Given an array $x$, this returns the $p$-norm
    /// $\|x\|\_p= (\sum_{j}|x_j|^{p})^{1/p}$.
    ///
    /// # Input
    /// ## p
    /// A positive number representing the exponent $p$ in the $p$-norm.
    /// ## array
    /// The vector $x$ of real numbers whose $p$-norm is to be returned.
    ///
    /// # Output
    /// The $p$-norm $\|x\|_p$.
    ///
    /// # Remarks
    /// This function defines a norm only when `p >= 1.0` or `Length(array)` is
    /// either 0 or 1. In the more general case, this function fails the
    /// triangle inequality.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Math.PNormalized
    function PNorm(p : Double, array : Double[]) : Double {
        if p <= 0.0 {
            fail $"PNorm failed. `p` must be a positive real number, but was {p}.";
        }

        mutable norm = 0.0;

        for element in array {
            set norm += PowD(AbsD(element), p);
        }

        return PowD(norm, 1.0 / p);
    }


    /// # Summary
    /// Returns the squared 2-norm of a vector.
    ///
    /// # Description
    /// Returns the squared 2-norm of a vector; that is, given an input
    /// $\vec{x}$, returns $\sum_i x_i^2$.
    ///
    /// # Input
    /// ## array
    /// The vector whose squared 2-norm is to be returned.
    ///
    /// # Output
    /// The squared 2-norm of `array`.
    function SquaredNorm(array : Double[]) : Double {
        mutable ret = 0.0;
        for element in array {
            set ret += element * element;
        }
        return ret;
    }


    /// # Summary
    /// Normalizes a vector of real numbers according to the p-norm for a given
    /// p.
    ///
    /// # Description
    /// That is, given an array $x$ of type `Double[]`, this returns an array where
    /// all elements are divided by the $p$-norm $\|x\|_p$.
    ///
    /// # Input
    /// ## p
    /// The exponent $p$ in the $p$-norm.
    /// ## array
    /// The vector $x$ to be normalized.
    ///
    /// # Output
    /// The array $x$ normalized by the $p$-norm $\|x\|_p$.
    ///
    /// # Remarks
    /// This function defines a norm only when `p >= 1.0` or `Length(array)` is
    /// either 0 or 1. In the more general case, this function fails the
    /// triangle inequality.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Math.PNorm
    function PNormalized(p : Double, array : Double[]) : Double[] {
        let nElements = Length(array);
        let norm = PNorm(p, array);

        if norm == 0.0 {
            return array;
        } else {
            mutable output = [0.0, size=nElements];

            for idx in 0 .. nElements - 1 {
                set output w/= idx <- array[idx] / norm;
            }

            return output;
        }
    }


    /// # Summary
    /// Returns the factorial of a given number.
    ///
    /// # Description
    /// Returns the factorial of a given nonnegative integer $n$, where $n \le 20$.
    ///
    /// # Input
    /// ## n
    /// The number to take the factorial of.
    ///
    /// # Output
    /// The factorial of `n`.
    ///
    /// # Remarks
    /// For inputs greater than 20, please use @"Microsoft.Quantum.Math.FactorialL".
    ///
    /// # See Also
    /// - Microsoft.Quantum.Math.FactorialL
    function FactorialI(n : Int) : Int {
        mutable an = 1;
        mutable x = 1;

        Fact(n >= 0, "The factorial is not defined for negative inputs.");
        Fact(n < 21, "The largest factorial that be stored as an Int is 20!. Use FactorialL or ApproximateFactorial.");

        if n == 0 {
            return x;
        } else {
            set an = n;
        }

        for i in  1 .. an {
            set x *= i;
        }

        return x;
    }


    /// # Summary
    /// Returns an approximate factorial of a given number.
    ///
    /// # Description
    /// Returns the factorial as `Double`, given an input of $n$ as a `Double`.
    /// The domain of inputs for this function is `AbsD(n) < 170.0`.
    ///
    /// # Remarks
    /// For $n \ge 10$, this function uses the Ramanujan approximation with a
    /// relative error to the order of $1 / n^5$.
    ///
    /// # Input
    /// ## n
    /// The number to take the approximate factorial of.
    ///
    /// # Output
    /// The approximate factorial of `n`.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Math.FactorialI
    /// - Microsoft.Quantum.Math.FactorialL
    function ApproximateFactorial(n : Int) : Double {
        Fact(n >= 0, "The factorial is not defined for negative inputs.");
        Fact(n < 170, "The largest approximate factorial that be stored as an Double is 169!. Use FactorialL.");

        // For small enough n, use the exact factorial instead.
        if n < 10 {
            return IntAsDouble(FactorialI(n));
        }

        let absN = IntAsDouble(n);

        let a = Sqrt(2.0 * PI() * absN);
        let b = (absN / E()) ^ absN;
        let c = E() ^ (1.0 / (12.0 * absN) - (1.0 / (360.0 * (absN ^ 3.0))));
        return a * b * c;
    }

    /// # Summary
    /// Returns the double factorial of a given integer.
    ///
    /// # Input
    /// ## n
    /// The number to take the double factorial of.
    ///
    /// # Output
    /// The double factorial of the provided input.
    ///
    /// # Remarks
    /// The double factorial $n!!$ of $n$ is defined as
    /// $n \times (n - 2) \times \cdots \times k$, where $k \in {1, 2}$. For example,
    /// $7!! = 7 \times 5 \times 3 \times 1$.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Math.ApproximateFactorial
    /// - Microsoft.Quantum.Math.FactorialI
    /// - Microsoft.Quantum.Math.FactorialL
    internal function DoubleFactorialL(n : Int) : BigInt {
        Fact(n >= 0, "The double factorial is not defined for negative inputs.");
        mutable acc = 1L;

        for i in (n % 2 == 0 ? 2 | 1)..2..AbsI(n) {
            set acc *= IntAsBigInt(i);
        }

        return acc;
    }


    /// # Summary
    /// Returns the factorial of a given integer.
    ///
    /// # Input
    /// ## n
    /// The number to take the factorial of.
    ///
    /// # Output
    /// The factorial of the provided input.
    ///
    /// # Remarks
    /// This function returns exact factorials for arbitrary-size integers,
    /// using a recursive decomposition into double-factorials ($n!!$).
    /// In particular, if $n = 2k + 1$ for $k \in \mathbb{N}$, then:
    /// $$
    ///     n! = n!! \times k! \times 2^k,
    /// $$
    /// where $k!$ can be computed recursively. If $n$ is even, then we can
    /// begin the recursion by computing $n! = n \times (n - 1)!$.
    ///
    ///
    /// # See Also
    /// - Microsoft.Quantum.Math.ApproximateFactorial
    /// - Microsoft.Quantum.Math.FactorialI
    function FactorialL(n : Int) : BigInt {
        if n < 0 {
            fail "The factorial is not defined for negative inputs.";
        }

        let absN = AbsI(n);

        if absN == 0 or absN == 1 {
            return 1L;
        }

        // If n is even, recurse on n - 1 so that we know we're starting at
        // an odd number.
        if n % 2 == 0 {
            return IntAsBigInt(n) * FactorialL(n - 1);
        }

        // At this point, we know that 𝑛 is an odd number >= 3.
        // Our approach will be to use that for 𝑛 = 2𝑘 + 1,
        // 𝑛! = 𝑛!! * 𝑘! * 2^𝑘.
        let k = absN / 2;
        return DoubleFactorialL(absN) * FactorialL(k) * (1L <<< k);
    }


    /// # Summary
    /// Returns the natural logarithm of the gamma function (aka the log-gamma
    /// function).
    ///
    /// # Description
    /// The gamma function $\Gamma(x)$ generalizes the factorial function
    /// to the positive real numbers and is defined as
    /// $$
    /// \begin{align}
    ///     \Gamma(x) \mathrel{:=} \int_0^{\infty} t^{x - 1} e^{-t} dt.
    /// \end{align}
    /// $$
    ///
    /// The gamma function has the property that for all positive real numbers
    /// $x$, $\Gamma(x + 1) = x \Gamma(x)$, such that the factorial function
    /// is a special case of $\Gamma$,
    /// $n! = \Gamma(n + 1)$ for all natural numbers $n$.
    ///
    /// # Input
    /// ## x
    /// The point $x$ at which the log-gamma function is to be evaluated.
    ///
    /// # Output
    /// The value $\ln \Gamma(x)$.
    function LogGammaD(x : Double) : Double {
        // Here, we use the approximation described in Numerical Recipes in C.
        let coefficients = [
            57.1562356658629235, -59.5979603554754912,
            14.1360979747417471, -0.491913816097620199, .339946499848118887e-4,
            .465236289270485756e-4, -.983744753048795646e-4, .158088703224912494e-3,
            -.210264441724104883e-3, .217439618115212643e-3, -.164318106536763890e-3,
            .844182239838527433e-4, -.261908384015814087e-4, .368991826595316234e-5
        ];

        Fact(x > 0.0, "Γ(x) not defined for x <= 0.");

        mutable y = x;
        let tmp = x + 5.2421875000000000;

        mutable acc = 0.99999999999999709;
        for coeff in coefficients {
            set y += 1.0;
            set acc += coeff / y;
        }


        return Log(2.506628274631000 * acc / x) + ((x + 0.5) * Log(tmp) - tmp);
    }

    /// # Summary
    /// Returns the approximate natural logarithm of the factorial of a given
    /// integer.
    ///
    /// # Input
    /// ## n
    /// The number to take the log-factorial of.
    ///
    /// # Output
    /// The natural logarithm of the factorial of the provided input.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Math.ApproximateFactorial
    /// - Microsoft.Quantum.Math.FactorialI
    /// - Microsoft.Quantum.Math.FactorialL
    function LogFactorialD(n : Int) : Double {
        return LogGammaD(IntAsDouble(n) + 1.0);
    }

    /// # Summary
    /// Returns the binomial coefficient of two integers.
    ///
    /// # Description
    /// Given two integers $n$ and $k$, returns the binomial coefficient
    /// $(n k)$, also known as $n$-choose-$k$.
    ///
    /// # Input
    /// ## n
    /// The first of the two integers to compute the binomial coefficient of.
    /// ## k
    /// The second of the two integers to compute the binomial coefficient of.
    ///
    /// # Output
    /// The binomial coefficient $(n k)$.
    function Binom(n : Int, k : Int) : Int {
        // Here, we use the approximation described in Numerical Recipes in C.
        if n < 171 {
            return Floor(0.5 + ApproximateFactorial(n) / (ApproximateFactorial(k) * ApproximateFactorial(n - k)));
        } else {
            return Floor(0.5 + ExpD(LogFactorialD(n) - LogFactorialD(k) - LogFactorialD(n - k)));
        }
    }

    /// # Summary
    /// Returns a binomial coefficient of the form "½-choose-k."
    ///
    /// # Description
    /// Given an integer $k$, returns the binomial coefficient
    /// $(\frac{1}{2} k)$, also known as $\frac{1}{2}$-choose-$k$.
    ///
    /// # Input
    /// ## k
    /// The integer to compute the half-integer binomial coefficient of.
    ///
    /// # Output
    /// The binomial coefficient $(\frac{1}{2} k)$.
    function HalfIntegerBinom(k : Int) : Double {
        let numerator = IntAsDouble(Binom(2 * k, k)) * IntAsDouble(k % 2 == 0 ? -1 | +1);
        return numerator / (2.0 ^ IntAsDouble(2 * k) * IntAsDouble(2 * k - 1));
    }

}
