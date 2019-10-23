// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {
    open Microsoft.Quantum.Math;

    /// # Deprecated
    /// This function has been removed.
    @Deprecated("Microsoft.Quantum.Math.AbsI")
    function IntAbs(input : Int) : Int {
        return AbsI(input);
    }

    /// # Deprecated
    /// This function has been removed.
    @Deprecated("Microsoft.Quantum.Math.MaxI")
    function IntMax(a : Int, b : Int) : Int {
        return MaxI(a, b);
    }

    /// # Deprecated
    /// Please use @"microsoft.quantum.math.pnormalized".
    @Deprecated("Microsoft.Quantum.Math.PNormalized")
    function PNormalize(p : Double, array : Double[]) : Double[] {
        return PNormalized(p, array);
    }

    /// # Deprecated
    /// Please use @"microsoft.quantum.math.extendedgreatestcommondivisori".
    @Deprecated("Microsoft.Quantum.Math.ExtendedGreatestCommonDivisorI")
    function ExtendedGCD (a : Int, b : Int) : (Int, Int) {
        return ExtendedGreatestCommonDivisorI(a, b);
    }

    /// # Deprecated
    /// Please use @"microsoft.quantum.math.greatestcommondivisori".
    @Deprecated("Microsoft.Quantum.Math.GreatestCommonDivisorI")
    function GCD (a : Int, b : Int) : Int {
        return GreatestCommonDivisorI(a, b);
    }

    /// # Deprecated
    /// Please use @"microsoft.quantum.math.greatestcommondivisori".
    @Deprecated("Microsoft.Quantum.Math.IsCoprimeI")
    function IsCoprime (a : Int, b : Int) : Bool {
        return IsCoprimeI(a, b);
    }

    /// # Deprecated
    /// Please use @"microsoft.quantum.math.inversemodi".
    @Deprecated("Microsoft.Quantum.Math.InverseModI")
    function InverseMod (a : Int, modulus : Int) : Int {
        return InverseModI(a, modulus);
    }

    /// # Deprecated
    /// Please use @"microsoft.quantum.math.modulusi".
    @Deprecated("Microsoft.Quantum.Math.ModulusI")
    function Modulus (value : Int, modulus : Int) : Int {
        return ModulusI(value, modulus);
    }

    /// # Deprecated
    /// Please use @"microsoft.quantum.math.expmodi".
    @Deprecated("Microsoft.Quantum.Math.ExpModI")
    function ExpMod (expBase : Int, power : Int, modulus : Int) : Int {
        return ExpModI(expBase, power, modulus);
    }

    /// # Deprecated
    /// Please use @"microsoft.quantum.math.continuedfractionconvergenti".
    @Deprecated("Microsoft.Quantum.Math.ContinuedFractionConvergentI")
    function ContinuedFractionConvergent (fraction : Fraction, denominatorBound : Int) : Fraction {
        return ContinuedFractionConvergentI(fraction, denominatorBound);
    }

    /// # Deprecated
    /// Please use @"microsoft.quantum.math.bitsizei".
    @Deprecated("Microsoft.Quantum.Math.BitSizeI")
    function BitSize(a : Int) : Int {
        return BitSizeI(a);
    }

}

namespace Microsoft.Quantum.Math {

    /// # Deprecated.
    /// Please use @"microsoft.quantum.math.complexpolarascomplex".
    @Deprecated("Microsoft.Quantum.Math.ComplexPolarAsComplex")
    function ComplexPolarToComplex(input : ComplexPolar) : Complex {
        return ComplexPolarAsComplex(input);
    }

    /// # Deprecated.
    /// Please use @"microsoft.quantum.math.complexascomplexpolar".
    @Deprecated("Microsoft.Quantum.Math.ComplexAsComplexPolar")
    function ComplexToComplexPolar(input : Complex) : ComplexPolar {
        return ComplexAsComplexPolar(input);
    }

}