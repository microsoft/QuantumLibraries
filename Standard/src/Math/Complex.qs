// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Math {

    /// # Summary
	/// Represents a complex number in polar form.
	///
    /// The polar representation of a complex number is $c=r e^{i t}$.
    ///
    /// # Input
    /// ## First Parameter
    /// `Double` is absolute value $r \ge 0$.
    /// ## Second Parameter
    /// `Double` is the phase $t \in \mathbb R$.
    newtype ComplexPolar = (Magnitude: Double, Argument: Double);

    /// # Summary
    /// Returns the squared absolute value of a complex number of type
    /// `Complex`.
    ///
    /// # Input
    /// ## input
    /// Complex number $c = x + i y$.
    ///
    /// # Output
    /// Squared absolute value $|c|^2 = x^2 + y^2$.
    function AbsSquaredComplex (input : Complex) : Double {
        let (real, imaginary) = input!;
        return real * real + imaginary * imaginary;
    }

    /// # Summary
    /// Returns the absolute value of a complex number of type
    /// `Complex`.
    ///
    /// # Input
    /// ## input
    /// Complex number $c = x + i y$.
    ///
    /// # Output
    /// Absolute value $|c| = \sqrt{x^2 + y^2}$.
    function AbsComplex (input : Complex) : Double {
        return Sqrt(AbsSquaredComplex(input));
    }

    /// # Summary
    /// Returns the phase of a complex number of type
    /// `Complex`.
    ///
    /// # Input
    /// ## input
    /// Complex number $c = x + i y$.
    ///
    /// # Output
    /// Phase $\text{Arg}[c] = \text{ArcTan}(y,x) \in (-\pi,\pi]$.
    function ArgComplex (input : Complex) : Double {
        let (real, imaginary) = input!;
        return ArcTan2(imaginary, real);
    }

    /// # Summary
    /// Returns the squared absolute value of a complex number of type
    /// `ComplexPolar`.
    ///
    /// # Input
    /// ## input
    /// Complex number $c = r e^{i t}$.
    ///
    /// # Output
    /// Squared absolute value $|c|^2 = r^2$.
    function AbsSquaredComplexPolar (input : ComplexPolar) : Double {
        let (abs, arg) = input!;
        return abs * abs;
    }
    
    
    /// # Summary
    /// Returns the absolute value of a complex number of type
    /// `ComplexPolar`.
    ///
    /// # Input
    /// ## input
    /// Complex number $c = r e^{i t}$.
    ///
    /// # Output
    /// Absolute value $|c| = r$.
    function AbsComplexPolar (input : ComplexPolar) : Double {
        return input::Magnitude;
    }
    
    
    /// # Summary
    /// Returns the phase of a complex number of type
    /// `ComplexPolar`.
    ///
    /// # Input
    /// ## input
    /// Complex number $c = r e^{i t}$.
    ///
    /// # Output
    /// Phase $\text{Arg}[c] = t$.
    function ArgComplexPolar (input : ComplexPolar) : Double {
        return input::Argument;
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
    function NegationC(input : Complex) : Complex {
        let (re, im) = input!;
        return Complex(-re, -im);
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
    function NegationCP(input : ComplexPolar) : ComplexPolar {
        return ComplexPolar(input::Magnitude, input::Argument + PI());
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
    function PlusC(a : Complex, b : Complex) : Complex {
        let ((reA, imA), (reB, imB)) = (a!, b!);
        return Complex(reA + reB, imA + imB);
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
    function PlusCP(a : ComplexPolar, b : ComplexPolar) : ComplexPolar {
        return ComplexAsComplexPolar(
            PlusC(
                ComplexPolarAsComplex(a),
                ComplexPolarAsComplex(b)
            )
        );
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
    function MinusC(a : Complex, b : Complex) : Complex {
        return PlusC(a, NegationC(b));
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
    function MinusCP(a : ComplexPolar, b : ComplexPolar) : ComplexPolar {
        return PlusCP(a, NegationCP(b));
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
    function TimesC(a : Complex, b : Complex) : Complex {
        let ((reA, imA), (reB, imB)) = (a!, b!);
        return Complex(
            reA * reB - imA * imB,
            reA * imB + imA * reB
        );
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
    function TimesCP(a : ComplexPolar, b : ComplexPolar) : ComplexPolar {
        return ComplexPolar(
            a::Magnitude * b::Magnitude,
            a::Argument + b::Argument
        );
    }

    /// # Summary
    /// Returns a number raised to a given power.
    ///
    /// # Input
    /// ## base
    /// The number $a$ that is to be raised.
    /// ## power
    /// The power $b$ to which $a$ should be raised.
    ///
    /// # Output
    /// The power $a^b$
    function PowC(a : Complex, power : Double) : Complex {
        return ComplexPolarAsComplex(
            PowCP(ComplexAsComplexPolar(a), power)
        );
    }

    /// # Summary
    /// Returns a number raised to a given power.
    ///
    /// # Input
    /// ## base
    /// The number $a$ that is to be raised.
    /// ## power
    /// The power $b$ to which $a$ should be raised.
    ///
    /// # Output
    /// The power $a^b$
    function PowCP(a : ComplexPolar, power : Double) : ComplexPolar {
        return ComplexPolar(
            PowD(a::Magnitude, power),
            a::Argument * power
        );
    }

}
