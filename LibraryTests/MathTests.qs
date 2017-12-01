// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Primitive;
    open Microsoft.Quantum.Extensions.Math;

    function NativeFnsAreCallableTest() : () {
        let arg = PI() / 2.0;
        AssertAlmostEqual(Sin(arg), 1.0);
        AssertAlmostEqual(Cos(arg), 0.0);

        let arcArg = 1.0;
        AssertAlmostEqual(ArcCos(arcArg), 0.0);
        AssertAlmostEqual(ArcSin(arcArg), arg);
    }

    function RealModTest() : () {
        AssertAlmostEqual(RealMod(5.5 * PI(), 2.0 * PI(), 0.0), 1.5 * PI());
        AssertAlmostEqual(RealMod(0.5 * PI(), 2.0 * PI(), -PI() / 2.0), 0.5 * PI());
    }

    function ArcHyperbolicFnsTest() : () {
        // These tests were generated using NumPy's implementations
        // of the inverse hyperbolic functions.
        AssertAlmostEqual(ArcTanh(0.3), 0.30951960420311175);
        AssertAlmostEqual(ArcCosh(1.3), 0.75643291085695963);
        AssertAlmostEqual(ArcSinh(-0.7), -0.65266656608235574);
    }

    function ExtendedGCDTestHelper(  a : Int , b : Int, gcd : Int ) : () {
        Message($"Testing {a}, {b}, {gcd} ");
        let (u,v) = ExtendedGCD(a,b);
        let expected = AbsI(gcd);
        let actual = AbsI(u*a+v*b);
        AssertIntEqual( expected, actual,
            $"Expected absolute value of gcd to be {expected}, got {actual}");
    }

    function ExtendedGCDTest() : ()
    {
        let testTuples = [ (1,1,1); (1,-1,1); (-1,1,1); (-1,-1,1); (5,7,1); (-5,7,1); (3,15,3) ];
        Ignore(Map(ExtendedGCDTestHelper, testTuples));
    }
}
