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

    @Test("QuantumSimulator")
    function NegationIsCorrect() : Unit {
        EqualityFactI(NegationI(42), -42, "NegationI returned wrong output.");
        NearEqualityFactD(NegationD(42.0), -42.0);
        EqualityFactL(NegationL(42L), -42L, "NegationI returned wrong output.");
        NearEqualityFactC(NegationC(Complex(-1.0, -2.0)), Complex(1.0, 2.0));
        NearEqualityFactCP(
            NegationCP(ComplexPolar(1.0, 5.0 * PI() / 4.0)),
            ComplexPolar(1.0, PI() / 4.0)
        );
    }

    @Test("QuantumSimulator")
    function PlusIsCorrect() : Unit {
        EqualityFactI(PlusI(-72, 32), -40, "PlusI returned wrong output.");
        NearEqualityFactD(PlusD(-72.0, 32.0), -40.0);
        EqualityFactL(PlusL(-72L, 32L), -40L, "PlusL returned wrong output.");
        NearEqualityFactC(PlusC(ONE_C(), ONE_C()), TWO_C());
        NearEqualityFactCP(PlusCP(ONE_CP(), ONE_CP()), TWO_CP());
    }

    @Test("QuantumSimulator")
    function MinusIsCorrect() : Unit {
        EqualityFactI(MinusI(72, 32), 40, "MinusI returned wrong output.");
        NearEqualityFactD(MinusD(72.0, 32.0), 40.0);
        EqualityFactL(MinusL(72L, 32L), 40L, "MinusL returned wrong output.");
        NearEqualityFactC(MinusC(TWO_C(), ONE_C()), ONE_C());
        NearEqualityFactCP(MinusCP(TWO_CP(), ONE_CP()), ONE_CP());
    }

    @Test("QuantumSimulator")
    function TimesIsCorrect() : Unit {
        EqualityFactI(TimesI(-10, -4), 40, "TimesI returned wrong output.");
        NearEqualityFactD(TimesD(-10.0, -4.0), 40.0);
        EqualityFactL(TimesL(-10L, -4L), 40L, "TimesL returned wrong output.");
        NearEqualityFactC(TimesC(TWO_C(), ONE_C()), TWO_C());
        NearEqualityFactCP(TimesCP(TWO_CP(), ONE_CP()), TWO_CP());
    }

    @Test("QuantumSimulator")
    function DividedByIsCorrect() : Unit {
        EqualityFactI(DividedByI(-40, -4), 10, "DividedByI returned wrong output.");
        NearEqualityFactD(DividedByD(-40.0, -4.0), 10.0);
        EqualityFactL(DividedByL(-40L, -4L), 10L, "DividedByL returned wrong output.");
        NearEqualityFactC(DividedByC(TWO_PI_4_C(), TWO_C()), PI_4_C());
        NearEqualityFactCP(DividedByCP(TWO_PI_4_CP(), TWO_CP()), PI_4_CP());
    }

    @Test("QuantumSimulator")
    function ModIsCorrect() : Unit {
        EqualityFactI(ModI(17, 5), 2, "ModI returned wrong output.");
        EqualityFactL(ModL(17L, 5L), 2L, "ModL returned wrong output.");
    }

    @Test("QuantumSimulator")
    function PowIsCorrect() : Unit {
        EqualityFactI(PowI(3, 8), 6561, "PowI returned wrong output.");
        NearEqualityFactD(PowD(2.0, 7.3), 157.58648490814928441592231285347);
        EqualityFactL(PowL(17L, 19), 239072435685151324847153L, "PowL returned wrong output.");
        NearEqualityFactC(PowC(PI_4_C(), TWO_C()), PI_2_C());
        NearEqualityFactCP(PowCP(PI_4_CP(), TWO_CP()), PI_2_CP());
    }

}
