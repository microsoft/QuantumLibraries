// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {
    open Microsoft.Quantum.Primitive;
    open Microsoft.Quantum.Extensions.Math;

	function Lg(input: Double) : Double
	{
		// Fully-qualified name is required because Log also appears in Primitives
		return Microsoft.Quantum.Extensions.Math.Log(input) / LogOf2();
	}

    // TODO: use generics for Double, Bool, Result
    /// summary:
    ///    Given an array of integers, returns the largest element.
    function Max(value : Int[]) : Int 
    {
        mutable max = value[0];
        let nTerms = Length(value);
        for(idx in 0..nTerms - 1)
        {
            if (value[idx]> max){
                set max = value[idx];
            }
        }
        return max;
    }

    /// summary:
    ///    Given an array of integers, returns the smallest element.
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
    
    /// summary:
    ///    Computes the modulus between two real numbers.
    /// example:
    ///     ```qflat
    ///         // Returns 3 p / 2.
    ///         let y = RealMod(5.5 * PI(), 2.0 * PI(), 0.0)
    ///         // Returns p / 2.
    ///         let z = RealMod(0.5 * PI(), 2.0 * PI(), -PI() / 2.0)
    ///     ```
    function RealMod(value: Double, modulo: Double, minValue: Double) : Double {
        let fractionalValue = 2.0 * PI() * ((value - minValue) / modulo - 0.5 );
        let cosFracValue = Cos(fractionalValue);
        let sinFracValue = Sin(fractionalValue);
        let moduloValue = 0.5 + ArcTan2(sinFracValue, cosFracValue) / ( 2.0 * PI() );
        let output = moduloValue * modulo + minValue;
        return output;
    }

    // NB: .NET's Math library does not provide hyperbolic arcfunctions.
	function ArcCosh(x : Double) : Double
	{
		// Fully-qualified name is required because Log also appears in Primitives
		return Microsoft.Quantum.Extensions.Math.Log(x + Sqrt(x * x - 1.0));
	}

	function ArcSinh(x : Double) : Double
	{
		// Fully-qualified name is required because Log also appears in Primitives
		return Microsoft.Quantum.Extensions.Math.Log(x + Sqrt(x * x + 1.0));
	}

	function ArcTanh(x : Double) : Double
	{
		// Fully-qualified name is required because Log also appears in Primitives
		return Microsoft.Quantum.Extensions.Math.Log((1.0 + x) / (1.0 - x)) * 0.5;
	}

}
