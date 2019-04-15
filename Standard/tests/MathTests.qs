// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.
namespace Microsoft.Quantum.Canon {
    open Microsoft.Quantum.Math;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Diagnostics;
    open Microsoft.Quantum.Arrays;

    function NativeFnsAreCallableTest () : Unit {
        
        let arg = PI() / 2.0;
        NearEqualityFact(Sin(arg), 1.0);
        NearEqualityFact(Cos(arg), 0.0);
        let arcArg = 1.0;
        NearEqualityFact(ArcCos(arcArg), 0.0);
        NearEqualityFact(ArcSin(arcArg), arg);
    }
    
    
    function RealModTest () : Unit {
        
        NearEqualityFact(RealMod(5.5 * PI(), 2.0 * PI(), 0.0), 1.5 * PI());
        NearEqualityFact(RealMod(0.5 * PI(), 2.0 * PI(), -PI() / 2.0), 0.5 * PI());
    }
    
    
    function ArcHyperbolicFnsTest () : Unit {
        
        // These tests were generated using NumPy's implementations
        // of the inverse hyperbolic functions.
        NearEqualityFact(ArcTanh(0.3), 0.30951960420311175);
        NearEqualityFact(ArcCosh(1.3), 0.75643291085695963);
        NearEqualityFact(ArcSinh(-0.7), -0.65266656608235574);
    }
    
    
    function ExtendedGCDTestHelper (a : Int, b : Int, gcd : Int) : Unit {
        
        Message($"Testing {a}, {b}, {gcd} ");
        let (u, v) = ExtendedGCD(a, b);
        let expected = AbsI(gcd);
        let actual = AbsI(u * a + v * b);
        EqualityFactI(expected, actual, $"Expected absolute value of gcd to be {expected}, got {actual}");
    }
    
    
    function ExtendedGCDTest () : Unit {
        
        let testTuples = [(1, 1, 1), (1, -1, 1), (-1, 1, 1), (-1, -1, 1), (5, 7, 1), (-5, 7, 1), (3, 15, 3)];
        Ignore(Mapped(ExtendedGCDTestHelper, testTuples));
    }
    
    
    function BitSizeTest () : Unit {
        
        EqualityFactI(BitSize(3), 2, $"BitSize(3) must be 2");
        EqualityFactI(BitSize(7), 3, $"BitSize(7) must be 2");
    }
    
    
    function ExpModTest () : Unit {
        
        // this test is generated using Mathematica PowerMod function
        let result = ExpMod(5, 4611686018427387903, 7);
        EqualityFactI(result, 6, $"The result must be 6, got {result}");
    }
    
    
    function ContinuedFractionConvergentTestHelper (numerator : Int, denominator : Int) : Unit {
        
        let bitSize = 2 * BitSize(denominator);
        let numeratorDyadic = (numerator * 2 ^ bitSize) / denominator;
        let (u, v) = (ContinuedFractionConvergent(Fraction(numeratorDyadic, 2 ^ bitSize), denominator))!;
        EqualityFactB(AbsI(u) == numerator and AbsI(v) == denominator, true, $"The result must be ±{numerator}/±{denominator} got {u}/{v}");
    }
    
    
    function ContinuedFractionConvergentEdgeCaseTestHelper (numerator : Int, denominator : Int, bound : Int) : Unit {
        
        let (num, denom) = (ContinuedFractionConvergent(Fraction(numerator, denominator), bound))!;
        EqualityFactB(AbsI(num) == numerator and AbsI(denom) == denominator, true, $"The result must be ±{numerator}/±{denominator} got {num}/{denom}");
    }
    
    
    function ContinuedFractionConvergentTest () : Unit {
        
        let testTuples = [(29, 47), (17, 37), (15, 67)];
        Ignore(Mapped(ContinuedFractionConvergentTestHelper, testTuples));
        let edgeCaseTestTuples = [(1, 4, 512), (3, 4, 512)];
        Ignore(Mapped(ContinuedFractionConvergentEdgeCaseTestHelper, edgeCaseTestTuples));
    }
    
    
    function ComplexMathTest () : Unit {
        
        mutable complexCases = [(0.123, 0.321), (0.123, -0.321), (-0.123, 0.321), (-0.123, -0.321)];
        
        for (idxCases in IndexRange(complexCases)) {
            let (complexRe, complexIm) = complexCases[idxCases];
            let complexAbs = Sqrt(complexRe * complexRe + complexIm * complexIm);
            let complexArg = ArcTan2(complexIm, complexRe);
            let complex = Complex(complexRe, complexIm);
            let complexPolar = ComplexPolar(complexAbs, complexArg);
            NearEqualityFact(AbsSquaredComplex(complex), complexAbs * complexAbs);
            NearEqualityFact(AbsComplex(complex), complexAbs);
            NearEqualityFact(ArgComplex(complex), complexArg);
            NearEqualityFact(AbsSquaredComplexPolar(complexPolar), complexAbs * complexAbs);
            NearEqualityFact(AbsComplexPolar(complexPolar), complexAbs);
            NearEqualityFact(ArgComplexPolar(complexPolar), complexArg);
            let (x, y) = (ComplexPolarToComplex(complexPolar))!;
            NearEqualityFact(x, complexRe);
            NearEqualityFact(y, complexIm);
            let (r, t) = (ComplexToComplexPolar(complex))!;
            NearEqualityFact(r, complexAbs);
            NearEqualityFact(t, complexArg);
        }
    }
    
    
    function PNormTest () : Unit {
        
        mutable testCases = [
            (1.0, [-0.1, 0.2, 0.3], 0.6), 
            (1.5, [0.1, -0.2, 0.3], 0.43346228721136096815),
            (2.0, [0.1, 0.2, -0.3], 0.37416573867739413856),
            (3.0, [0.0, 0.0, -0.0], 0.0)
        ];
        
        for (idxTest in IndexRange(testCases)) {
            let (p, array, pNormExpected) = testCases[idxTest];
            NearEqualityFact(PNorm(p, array), pNormExpected);
            
            // if PNorm fails, PNormalized will fail.
            let arrayNormalized = PNormalized(p, array);
            
            for (idxCoeff in IndexRange(array)) {
                NearEqualityFact(array[idxCoeff] / pNormExpected, arrayNormalized[idxCoeff]);
            }
        }
    }
    
}


