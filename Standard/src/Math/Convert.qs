// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Math {

    /// # Summary
    /// Converts a real floating-point number to a complex number in its polar
    /// representation.
    ///
    /// # Input
    /// ## input
    /// The real number to be represented as a complex number.
    ///
    /// # Output
    /// A complex number representing the given input in terms of polar
    /// coordinates.
    function DoubleAsComplexPolar(input : Double) : ComplexPolar {
        return ComplexPolar(
            AbsD(input),
            input < 0.0 ? PI() | 0.0
        );
    }

    /// # Summary
    /// Converts a complex number of type `ComplexPolar` to a complex
    /// number of type `Complex`.
    ///
    /// # Input
    /// ## input
    /// Complex number $c = r e^{i t}$.
    ///
    /// # Output
    /// Complex number $c = x + i y$.
    function ComplexPolarAsComplex (input : ComplexPolar) : Complex {
        return Complex(
            input::Magnitude * Cos(input::Argument),
            input::Magnitude * Sin(input::Argument)
        );
    }

    /// # Summary
    /// Converts a complex number of type `Complex` to a complex
    /// number of type `ComplexPolar`.
    ///
    /// # Input
    /// ## input
    /// Complex number $c = x + i y$.
    ///
    /// # Output
    /// Complex number $c = r e^{i t}$.
    function ComplexAsComplexPolar (input : Complex) : ComplexPolar {
        return ComplexPolar(AbsComplex(input), ArgComplex(input));
    }

}
