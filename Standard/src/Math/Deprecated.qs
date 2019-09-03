// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {
    open Microsoft.Quantum.Math;
    open Microsoft.Quantum.Warnings;

    /// # Deprecated
    /// This function has been removed.
    function IntAbs(input : Int) : Int {
        _Removed("Microsoft.Quantum.Canon.IntAbs", "Microsoft.Quantum.Math.AbsI");
        return AbsI(input);
    }

    /// # Deprecated
    /// This function has been removed.
    function IntMax(a : Int, b : Int) : Int {
        _Removed("Microsoft.Quantum.Canon.IntMax", "Microsoft.Quantum.Math.MaxI");
        return MaxI(a, b);
    }

    /// # Deprecated
    /// Please use @"microsoft.quantum.math.pnormalized".
    function PNormalize(p : Double, array : Double[]) : Double[] {
        _Renamed("Microsoft.Quantum.Canon.PNormalize", "Microsoft.Quantum.Math.PNormalized");
        return PNormalized(p, array);
    }

    /// # Deprecated
    /// Please use @"microsoft.quantum.math.extendedgreatestcommondivisori".
    function ExtendedGCD (a : Int, b : Int) : (Int, Int) {
        _Renamed("Microsoft.Quantum.Canon.ExtendedGCD", "Microsoft.Quantum.Math.ExtendedGreatestCommonDivisorI");
        return ExtendedGreatestCommonDivisorI(a, b);
    }

    /// # Deprecated
    /// Please use @"microsoft.quantum.math.greatestcommondivisori".
    function GCD (a : Int, b : Int) : Int {
        _Renamed("Microsoft.Quantum.Canon.GCD", "Microsoft.Quantum.Math.GreatestCommonDivisorI");
        return GreatestCommonDivisorI(a, b);
    }

    /// # Deprecated
    /// Please use @"microsoft.quantum.math.greatestcommondivisori".
    function IsCoprime (a : Int, b : Int) : Bool {
        _Renamed("Microsoft.Quantum.Canon.IsCoprime", "Microsoft.Quantum.Math.IsCoprimeI");
        return IsCoprimeI(a, b);
    }

    /// # Deprecated
    /// Please use @"microsoft.quantum.math.inversemodi".
    function InverseMod (a : Int, modulus : Int) : Int {
        _Renamed("Microsoft.Quantum.Canon.InverseMod", "Microsoft.Quantum.Math.InverseModI");
        return InverseModI(a, modulus);
    }

    /// # Deprecated
    /// Please use @"microsoft.quantum.math.modulusi".
    function Modulus (value : Int, modulus : Int) : Int {
        _Renamed("Microsoft.Quantum.Canon.Modulus", "Microsoft.Quantum.Math.ModulusI");
        return ModulusI(value, modulus);
    }

    /// # Deprecated
    /// Please use @"microsoft.quantum.math.expmodi".
    function ExpMod (expBase : Int, power : Int, modulus : Int) : Int {
        _Renamed("Microsoft.Quantum.Canon.ExpMod", "Microsoft.Quantum.Math.ExpModI");
        return ExpModI(expBase, power, modulus);
    }

    /// # Deprecated
    /// Please use @"microsoft.quantum.math.continuedfractionconvergenti".
    function ContinuedFractionConvergent (fraction : Fraction, denominatorBound : Int) : Fraction {
        _Renamed("Microsoft.Quantum.Canon.ContinuedFractionConvergent", "Microsoft.Quantum.Math.ContinuedFractionConvergentI");
        return ContinuedFractionConvergent(fraction, denominatorBound);
    }

    /// # Deprecated
    /// Please use @"microsoft.quantum.math.bitsizei".
    function BitSize(a : Int) : Int {
        _Renamed("Microsoft.Quantum.Canon.BitSize", "Microsoft.Quantum.Math.BitSizeI");
        return BitSizeI(a);
    }

}

namespace Microsoft.Quantum.Math {
    open Microsoft.Quantum.Warnings;

    /// # Deprecated.
    /// Please use @"microsoft.quantum.math.complexpolarascomplex".
    function ComplexPolarToComplex(input : ComplexPolar) : Complex {
        _Renamed("Microsoft.Quantum.Math.ComplexPolarToComplex", "Microsoft.Quantum.Math.ComplexPolarAsComplex");
        return ComplexPolarAsComplex(input);
    }

    /// # Deprecated.
    /// Please use @"microsoft.quantum.math.complexascomplexpolar".
    function ComplexToComplexPolar(input : Complex) : ComplexPolar {
        _Renamed("Microsoft.Quantum.Math.ComplexToComplexPolar", "Microsoft.Quantum.Math.ComplexAsComplexPolar");
        return ComplexAsComplexPolar(input);
    }

}