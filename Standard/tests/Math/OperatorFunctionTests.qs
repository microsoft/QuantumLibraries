// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Tests {
    open Microsoft.Quantum.Diagnostics;
    open Microsoft.Quantum.Math;

    // Define some useful constants for complex arithmetic.
    function ONE_C() : Complex { return Complex(1.0, 0.0); }
    function ONE_CP() : ComplexPolar { return ComplexPolar(1.0, 0.0); }
    function TWO_C() : Complex { return Complex(2.0, 0.0); }
    function TWO_CP() : ComplexPolar { return ComplexPolar(2.0, 0.0); }
    function PI_2_C() : Complex { return ComplexPolarAsComplex(PI_2_CP()); }
    function PI_2_CP() : ComplexPolar { return ComplexPolar(1.0, PI() / 2.0); }
    function PI_4_C() : Complex { return ComplexPolarAsComplex(PI_4_CP()); }
    function PI_4_CP() : ComplexPolar { return ComplexPolar(1.0, PI() / 4.0); }
    function TWO_PI_4_C() : Complex { return ComplexPolarAsComplex(TWO_PI_4_CP()); }
    function TWO_PI_4_CP() : ComplexPolar { return ComplexPolar(2.0, PI() / 4.0); }

    // FIXME: expected and actual are flipped uniformly in this file.
    //        this has no effect other than making unit test failures slightly
    //        harder to read.

    @Test("QuantumSimulator")
    function NegationIsCorrect() : Unit {
        EqualityFactI(-42, NegationI(42), "NegationI returned wrong output.");
        NearEqualityFactD(-42.0, NegationD(42.0));
        EqualityFactL(-42L, NegationL(42L), "NegationI returned wrong output.");
        NearEqualityFactC(Complex(1.0, 2.0), NegationC(Complex(-1.0, -2.0)));
        NearEqualityFactCP(
            ComplexPolar(1.0, PI() / 4.0),
            NegationCP(ComplexPolar(1.0, 5.0 * PI() / 4.0))
        );
    }

    @Test("QuantumSimulator")
    function PlusIsCorrect() : Unit {
        EqualityFactI(-40, PlusI(-72, 32), "PlusI returned wrong output.");
        NearEqualityFactD(-40.0, PlusD(-72.0, 32.0));
        EqualityFactL(-40L, PlusL(-72L, 32L), "PlusL returned wrong output.");
        NearEqualityFactC(TWO_C(), PlusC(ONE_C(), ONE_C()));
        NearEqualityFactCP(TWO_CP(), PlusCP(ONE_CP(), ONE_CP()));
    }

    @Test("QuantumSimulator")
    function MinusIsCorrect() : Unit {
        EqualityFactI(40, MinusI(72, 32), "MinusI returned wrong output.");
        NearEqualityFactD(40.0, MinusD(72.0, 32.0));
        EqualityFactL(40L, MinusL(72L, 32L), "MinusL returned wrong output.");
        NearEqualityFactC(ONE_C(), MinusC(TWO_C(), ONE_C()));
        NearEqualityFactCP(ONE_CP(), MinusCP(TWO_CP(), ONE_CP()));
    }

    @Test("QuantumSimulator")
    function TimesIsCorrect() : Unit {
        EqualityFactI(40, TimesI(-10, -4), "TimesI returned wrong output.");
        NearEqualityFactD(40.0, TimesD(-10.0, -4.0));
        EqualityFactL(40L, TimesL(-10L, -4L), "TimesL returned wrong output.");
        NearEqualityFactC(TWO_C(), TimesC(TWO_C(), ONE_C()));
        NearEqualityFactCP(TWO_CP(), TimesCP(TWO_CP(), ONE_CP()));
    }

    @Test("QuantumSimulator")
    function DividedByIsCorrect() : Unit {
        EqualityFactI(10, DividedByI(-40, -4), "DividedByI returned wrong output.");
        NearEqualityFactD(10.0, DividedByD(-40.0, -4.0));
        EqualityFactL(10L, DividedByL(-40L, -4L), "DividedByL returned wrong output.");
        NearEqualityFactC(PI_4_C(), DividedByC(TWO_PI_4_C(), TWO_C()));
        NearEqualityFactCP(PI_4_CP(), DividedByCP(TWO_PI_4_CP(), TWO_CP()));
    }

    @Test("QuantumSimulator")
    function ModIsCorrect() : Unit {
        EqualityFactI(2, ModI(17, 5), "ModI returned wrong output.");
        EqualityFactL(2L, ModL(17L, 5L), "ModL returned wrong output.");
    }

    @Test("QuantumSimulator")
    function PowIsCorrect() : Unit {
        EqualityFactI(6561, PowI(3, 8), "PowI returned wrong output.");
        NearEqualityFactD(157.58648490814928441592231285347, PowD(2.0, 7.3));
        EqualityFactL(239072435685151324847153L, PowL(17L, 19), "PowL returned wrong output.");
        NearEqualityFactC(PI_2_C(), PowC(PI_4_C(), TWO_C()));
        NearEqualityFactCP(PI_2_CP(), PowCP(PI_4_CP(), TWO_CP()));
    }

}
