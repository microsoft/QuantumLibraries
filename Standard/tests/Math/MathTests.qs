// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.
namespace Microsoft.Quantum.Canon {
    open Microsoft.Quantum.Math;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Diagnostics;
    open Microsoft.Quantum.Arrays;

    @Test("QuantumSimulator")
    function NativeFnsAreCallable () : Unit {
        
        let arg = PI() / 2.0;
        NearEqualityFactD(Sin(arg), 1.0);
        NearEqualityFactD(Cos(arg), 0.0);
        let arcArg = 1.0;
        NearEqualityFactD(ArcCos(arcArg), 0.0);
        NearEqualityFactD(ArcSin(arcArg), arg);
    }
    
    @Test("QuantumSimulator")
    function RealModIsCorrect () : Unit {
        
        NearEqualityFactD(RealMod(5.5 * PI(), 2.0 * PI(), 0.0), 1.5 * PI());
        NearEqualityFactD(RealMod(0.5 * PI(), 2.0 * PI(), -PI() / 2.0), 0.5 * PI());
    }
    
    @Test("QuantumSimulator")
    function ArcHyperbolicFnsAreCorrect () : Unit {
        
        // These tests were generated using NumPy's implementations
        // of the inverse hyperbolic functions.
        NearEqualityFactD(ArcTanh(0.3), 0.30951960420311175);
        NearEqualityFactD(ArcCosh(1.3), 0.75643291085695963);
        NearEqualityFactD(ArcSinh(-0.7), -0.65266656608235574);
    }
    
    function ExtendedGreatestCommonDivisorITestHelper (a : Int, b : Int, gcd : Int) : Unit {
        
        Message($"Testing {a}, {b}, {gcd} ");
        let (u, v) = ExtendedGreatestCommonDivisorI(a, b);
        let expected = AbsI(gcd);
        let actual = AbsI(u * a + v * b);
        EqualityFactI(expected, actual, $"Expected absolute value of gcd to be {expected}, got {actual}");
    }
    
    @Test("QuantumSimulator")
    function ExtendedGreatestCommonDivisorIIsCorrect () : Unit {
        
        let testTuples = [(1, 1, 1), (1, -1, 1), (-1, 1, 1), (-1, -1, 1), (5, 7, 1), (-5, 7, 1), (3, 15, 3)];
        Ignore(Mapped(ExtendedGreatestCommonDivisorITestHelper, testTuples));
    }

    @Test("QuantumSimulator")
    function GreatestCommonDivisorLIsCorrect() : Unit {
        EqualityFactL(
            GreatestCommonDivisorL(
                44958225298979240833230460209285719018635426448048959524915L,
                935899140510257015115572178639093206049354855486377228520740L
            ),
            5L,
            "GCD returned wrong result for BigInt inputs."
        );
    }
    
    @Test("QuantumSimulator")
    function BitSizeIsCorrect () : Unit {
        
        EqualityFactI(BitSizeI(3), 2, $"BitSizeI(3) must be 2");
        EqualityFactI(BitSizeI(7), 3, $"BitSizeI(7) must be 2");
    }

    @Test("QuantumSimulator")
    function CanComputeBitSizeFromLargeNumbers () : Unit {
        for (k in 1 .. 100) {
            let exp = 128 * k;
            EqualityFactI(BitSizeL(1L <<< exp), exp + 1, $"unexpected bitsize for exponent {exp} (k = {k})");
        }
    }
    
    @Test("QuantumSimulator")
    function ExpModIsCorrect () : Unit {
        
        // this test is generated using Mathematica PowerMod function
        let result = ExpModI(5, 4611686018427387903, 7);
        EqualityFactI(result, 6, $"The result must be 6, got {result}");
    }
     
    function ContinuedFractionConvergentTestHelper (numerator : Int, denominator : Int) : Unit {
        
        let bitSize = 2 * BitSizeI(denominator);
        let numeratorDyadic = (numerator * 2 ^ bitSize) / denominator;
        let (u, v) = (ContinuedFractionConvergentI(Fraction(numeratorDyadic, 2 ^ bitSize), denominator))!;
        EqualityFactB(AbsI(u) == numerator and AbsI(v) == denominator, true, $"The result must be ±{numerator}/±{denominator} got {u}/{v}");
    }
    
    function ContinuedFractionConvergentEdgeCaseTestHelper (numerator : Int, denominator : Int, bound : Int) : Unit {
        
        let (num, denom) = (ContinuedFractionConvergentI(Fraction(numerator, denominator), bound))!;
        EqualityFactB(AbsI(num) == numerator and AbsI(denom) == denominator, true, $"The result must be ±{numerator}/±{denominator} got {num}/{denom}");
    }
    
    @Test("QuantumSimulator")
    function ContinuedFractionConvergentIsCorrect () : Unit {
        
        let testTuples = [(29, 47), (17, 37), (15, 67)];
        Ignore(Mapped(ContinuedFractionConvergentTestHelper, testTuples));
        let edgeCaseTestTuples = [(1, 4, 512), (3, 4, 512)];
        Ignore(Mapped(ContinuedFractionConvergentEdgeCaseTestHelper, edgeCaseTestTuples));
    }
    
    @Test("QuantumSimulator")
    function ComplexMathIsCorrect () : Unit {
        
        mutable complexCases = [(0.123, 0.321), (0.123, -0.321), (-0.123, 0.321), (-0.123, -0.321)];
        
        for (idxCases in IndexRange(complexCases)) {
            let (complexRe, complexIm) = complexCases[idxCases];
            let complexAbs = Sqrt(complexRe * complexRe + complexIm * complexIm);
            let complexArg = ArcTan2(complexIm, complexRe);
            let complex = Complex(complexRe, complexIm);
            let complexPolar = ComplexPolar(complexAbs, complexArg);
            NearEqualityFactD(AbsSquaredComplex(complex), complexAbs * complexAbs);
            NearEqualityFactD(AbsComplex(complex), complexAbs);
            NearEqualityFactD(ArgComplex(complex), complexArg);
            NearEqualityFactD(AbsSquaredComplexPolar(complexPolar), complexAbs * complexAbs);
            NearEqualityFactD(AbsComplexPolar(complexPolar), complexAbs);
            NearEqualityFactD(ArgComplexPolar(complexPolar), complexArg);
            let (x, y) = (ComplexPolarAsComplex(complexPolar))!;
            NearEqualityFactD(x, complexRe);
            NearEqualityFactD(y, complexIm);
            let (r, t) = (ComplexAsComplexPolar(complex))!;
            NearEqualityFactD(r, complexAbs);
            NearEqualityFactD(t, complexArg);
        }
    }
    
    @Test("QuantumSimulator")
    function PNormIsCorrect () : Unit {
        
        mutable testCases = [
            (1.0, [-0.1, 0.2, 0.3], 0.6), 
            (1.5, [0.1, -0.2, 0.3], 0.43346228721136096815),
            (2.0, [0.1, 0.2, -0.3], 0.37416573867739413856),
            (3.0, [0.0, 0.0, -0.0], 0.0)
        ];
        
        for (idxTest in IndexRange(testCases)) {
            let (p, array, pNormExpected) = testCases[idxTest];
            NearEqualityFactD(PNorm(p, array), pNormExpected);
            
            // if PNorm fails, PNormalized will fail.
            let arrayNormalized = PNormalized(p, array);
            
            for (idxCoeff in IndexRange(array)) {
                NearEqualityFactD(array[idxCoeff] / pNormExpected, arrayNormalized[idxCoeff]);
            }
        }
    }

    @Test("QuantumSimulator")
    function SquaredNormIsCorrect() : Unit {
        NearEqualityFactD(SquaredNorm([2.0]), 4.0);
        NearEqualityFactD(SquaredNorm([1.0, 1.0]), 2.0);
        NearEqualityFactD(SquaredNorm([3.0, 4.0]), 25.0);
    }
}
