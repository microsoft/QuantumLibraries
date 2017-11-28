// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {
    open Microsoft.Quantum.Primitive;
    open Microsoft.Quantum.Extensions.Math;

    /// # Summary
    /// Computes the base-2 logarithm of a number.
    ///
    /// # Input
    /// ## input
    /// A real number $x$.
    ///
    /// # Output
    /// The base-2 logarithm $y = \log_2(x)$ such that $x = 2^y$.
	function Lg(input: Double) : Double
	{
		// Fully-qualified name is required because Log also appears in Primitives
		return Microsoft.Quantum.Extensions.Math.Log(input) / LogOf2();
	}

    // TODO: rewrite in terms of Fold(PairwiseMax, INTEGER_MIN(), _).
    /// # Summary
    /// Given an array of integers, returns the largest element.
    ///
    /// # Input
    /// ## values
    /// An array to take the maximum of.
    ///
    /// # Output
    /// The largest element of `values`.
    function Max(values : Int[]) : Int 
    {
        mutable max = values[0];
        let nTerms = Length(values);
        for(idx in 0..nTerms - 1)
        {
            if (values[idx]> max){
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
    function Min(value : Int[]) : Int 
    {
        mutable min = value[0];
        let nTerms = Length(value);
        for(idx in 0..nTerms - 1)
        {
            if (value[idx] < min) {
                set min = value[idx];
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
    /// # Example
    /// ```qsharp
    ///     // Returns 3 π / 2.
    ///     let y = RealMod(5.5 * PI(), 2.0 * PI(), 0.0)
    ///     // Returns -1.2, since +3.6 and -1.2 are 4.8 apart on the real line,
    ///     // which is a multiple of 2.4.
    ///     let z = RealMod(3.6, 2.4, -1.2);
    /// ```
    function RealMod(value: Double, modulo: Double, minValue: Double) : Double {
        let fractionalValue = 2.0 * PI() * ((value - minValue) / modulo - 0.5 );
        let cosFracValue = Cos(fractionalValue);
        let sinFracValue = Sin(fractionalValue);
        let moduloValue = 0.5 + ArcTan2(sinFracValue, cosFracValue) / ( 2.0 * PI() );
        let output = moduloValue * modulo + minValue;
        return output;
    }

    // NB: .NET's Math library does not provide hyperbolic arcfunctions.

    /// # Summary
    /// Computes the inverse hyperbolic cosine of a number.
    ///
    /// # Input
    /// ## x
    /// A real number $x$.
    ///
    /// # Output
    /// A real number $y$ such that $x = \cosh(y)$.
	function ArcCosh(x : Double) : Double
	{
		// Fully-qualified name is required because Log also appears in Primitives
		return Microsoft.Quantum.Extensions.Math.Log(x + Sqrt(x * x - 1.0));
	}

	/// # Summary
    /// Computes the inverse hyperbolic secant of a number.
    ///
    /// # Input
    /// ## x
    /// A real number $x$.
    ///
    /// # Output
    /// A real number $y$ such that $x = \sech(y)$.
	function ArcSinh(x : Double) : Double
	{
		// Fully-qualified name is required because Log also appears in Primitives
		return Microsoft.Quantum.Extensions.Math.Log(x + Sqrt(x * x + 1.0));
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
	function ArcTanh(x : Double) : Double
	{
		// Fully-qualified name is required because Log also appears in Primitives
		return Microsoft.Quantum.Extensions.Math.Log((1.0 + x) / (1.0 - x)) * 0.5;
	}

}
