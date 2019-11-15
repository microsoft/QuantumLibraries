// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Tests {
    open Microsoft.Quantum.Logical;
    open Microsoft.Quantum.Diagnostics;
    open Microsoft.Quantum.Math;

    function EqualTest() : Unit {
        Fact(EqualI(42, 42), "EqualI returned wrong output.");
        Fact(not EqualI(42, 73), "EqualI returned wrong output.");

        Fact(EqualL(42L, 42L), "EqualL returned wrong output.");
        Fact(not EqualL(42L, 730L), "EqualL returned wrong output.");

        Fact(EqualD(42.0, 42.0), "EqualD returned wrong output.");
        Fact(not EqualD(42.0, 730.0), "EqualD returned wrong output.");

        Fact(EqualR(One, One), "EqualR returned wrong output.");
        Fact(not EqualR(Zero, One), "EqualR returned wrong output.");

        Fact(EqualB(true, true), "EqualB returned wrong output.");
        Fact(not EqualB(true, false), "EqualB returned wrong output.");

        Fact(EqualC(Complex(1.0, 2.0), Complex(1.0, 2.0)), "EqualC returned wrong output.");
        Fact(not EqualC(Complex(1.0, 2.0), Complex(1.0, 73.0)), "EqualC returned wrong output.");

        Fact(EqualCP(ComplexPolar(1.0, 2.0), ComplexPolar(1.0, 2.0)), "EqualCP returned wrong output.");
        Fact(not EqualCP(ComplexPolar(1.0, 2.0), ComplexPolar(1.0, 73.0)), "EqualCP returned wrong output.");
    }

    function NotEqualTest() : Unit {
        Fact(not NotEqualI(42, 42), "NotEqualI returned wrong output.");
        Fact(NotEqualI(42, 73), "NotEqualI returned wrong output.");

        Fact(not NotEqualL(42L, 42L), "NotEqualL returned wrong output.");
        Fact(NotEqualL(42L, 730L), "NotEqualL returned wrong output.");

        Fact(not NotEqualD(42.0, 42.0), "NotEqualD returned wrong output.");
        Fact(NotEqualD(42.0, 730.0), "NotEqualD returned wrong output.");

        Fact(not NotEqualR(One, One), "NotEqualR returned wrong output.");
        Fact(NotEqualR(Zero, One), "NotEqualR returned wrong output.");

        Fact(not NotEqualB(true, true), "NotEqualB returned wrong output.");
        Fact(NotEqualB(true, false), "NotEqualB returned wrong output.");

        Fact(not NotEqualC(Complex(1.0, 2.0), Complex(1.0, 2.0)), "NotEqualC returned wrong output.");
        Fact(NotEqualC(Complex(1.0, 2.0), Complex(1.0, 73.0)), "NotEqualC returned wrong output.");

        Fact(not NotEqualCP(ComplexPolar(1.0, 2.0), ComplexPolar(1.0, 2.0)), "NotEqualCP returned wrong output.");
        Fact(NotEqualCP(ComplexPolar(1.0, 2.0), ComplexPolar(1.0, 73.0)), "NotEqualCP returned wrong output.");
    }

    function GreaterThanTest() : Unit {
        Fact(GreaterThanI(75, 32), "GreaterThanI returned wrong output.");
        Fact(not GreaterThanI(-13, 32), "GreaterThanI returned wrong output.");

        Fact(GreaterThanD(75.0, 32.0), "GreaterThanD returned wrong output.");
        Fact(not GreaterThanD(-13.0, 32.0), "GreaterThanD returned wrong output.");

        Fact(GreaterThanL(75L, 32L), "GreaterThanL returned wrong output.");
        Fact(not GreaterThanL(-13L, 32L), "GreaterThanL returned wrong output.");
    }

    function LessThanTest() : Unit {
        Fact(not LessThanI(75, 32), "LessThanI returned wrong output.");
        Fact(LessThanI(-13, 32), "LessThanI returned wrong output.");

        Fact(not LessThanD(75.0, 32.0), "LessThanD returned wrong output.");
        Fact(LessThanD(-13.0, 32.0), "LessThanD returned wrong output.");

        Fact(not LessThanL(75L, 32L), "LessThanL returned wrong output.");
        Fact(LessThanL(-13L, 32L), "LessThanL returned wrong output.");
    }

    function GreaterThanOrEqualTest() : Unit {
        Fact(GreaterThanOrEqualI(75, 75), "GreaterThanOrEqualI returned wrong output.");
        Fact(not GreaterThanOrEqualI(-13, 32), "GreaterThanOrEqualI returned wrong output.");

        Fact(GreaterThanOrEqualD(75.0, 75.0), "GreaterThanOrEqualD returned wrong output.");
        Fact(not GreaterThanOrEqualD(-13.0, 32.0), "GreaterThanOrEqualD returned wrong output.");

        Fact(GreaterThanOrEqualL(75L, 75L), "GreaterThanOrEqualL returned wrong output.");
        Fact(not GreaterThanOrEqualL(-13L, 32L), "GreaterThanOrEqualL returned wrong output.");
    }

    function LessThanOrEqualTest() : Unit {
        Fact(LessThanOrEqualI(75, 75), "LessThanOrEqualI returned wrong output.");
        Fact(not LessThanOrEqualI(32, -13), "LessThanOrEqualI returned wrong output.");

        Fact(LessThanOrEqualD(75.0, 75.0), "LessThanOrEqualD returned wrong output.");
        Fact(not LessThanOrEqualD(32.0, -13.0), "LessThanOrEqualD returned wrong output.");

        Fact(LessThanOrEqualL(75L, 75L), "LessThanOrEqualL returned wrong output.");
        Fact(not LessThanOrEqualL(32L, -13L), "LessThanOrEqualL returned wrong output.");
    }

    function NearlyEqualDTest() : Unit {
        Fact(NearlyEqualD(1.0, 1.0), "Exactly equal numbers marked as not nearly equal.");
        Fact(NearlyEqualDTest(1.0, 1.0 + 1e-15), "Nearly equal numbers marked as not nearly equal.");
        Fact(not NearlyEqualD(1.0, 1000.0), "Not nearly equal numbers marked as nearly equal.");
    }

}
