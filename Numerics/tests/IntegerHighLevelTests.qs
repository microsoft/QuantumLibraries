// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Tests {
    open Microsoft.Quantum.Diagnostics;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Arithmetic;
    open Microsoft.Quantum.Extensions.Math;

    @Test("ToffoliSimulator")
    operation MultiplyIExhaustiveTest() : Unit {
        ExhaustiveTestHelper2NonRegularArgs(IntegerMultiplicationRun(false, _, _, _, _, _));
    }

    @Test("ToffoliSimulator")
    operation SquareIExhaustiveTest() : Unit {
        ExhaustiveTestHelper1Arg(IntegerSquareRun(false, _, _, _));
    }

    @Test("ToffoliSimulator")
    operation DivideIExhaustiveTest() : Unit {
        ExhaustiveTestHelper2RegularArgs(IntegerDivisionRun);
    }

    @Test("ToffoliSimulator")
    operation SquareSIExhaustiveTest() : Unit {
        ExhaustiveTestHelper1Arg(IntegerSquareRun(true, _, _, _));
    }

    @Test("ToffoliSimulator")
    operation CompareGTSIExhaustiveTest() : Unit {
        ExhaustiveTestHelper2RegularArgs(IntegerGreaterThanRun(true, _, _, _, _));
    }

    @Test("ToffoliSimulator")
    operation MultiplySIExhaustiveTest() : Unit {
        ExhaustiveTestHelper2NonRegularArgs(IntegerMultiplicationRun(true, _, _, _, _, _));
    }

    @Test("ToffoliSimulator")
    operation ComputeReciprocalIExhaustiveTest() : Unit {
        ExhaustiveTestHelper1Arg(IntegerReciprocalRun(false, _, _, _));
    }

    internal operation IntegerGreaterThanRun(signed: Bool, a: Int, b: Int,
                                             n: Int, numCtrl: Int) : Unit {
        use aqs = Qubit[n];
        use bqs = Qubit[n];
        use result = Qubit();
        use ctrlqs = Qubit[numCtrl];
        ApplyXorInPlace(a, LittleEndian(aqs));
        ApplyXorInPlace(b, LittleEndian(bqs));
        if signed {
            CompareGTSI(
                SignedLittleEndian(LittleEndian(aqs)),
                SignedLittleEndian(LittleEndian(bqs)),
                result);
        } else {
            CompareGTI(LittleEndian(aqs),
                LittleEndian(bqs),
                result);
        }
        mutable asigned = a;
        mutable bsigned = b;
        if signed and a >= 2^(n-1) {
            set asigned = -2^n+a;
        }
        if signed and b >= 2^(n-1) {
            set bsigned = -2^n+b;
        }
        mutable res = asigned > bsigned;
        mutable resMeasured = M(result);
        EqualityFactB(res, resMeasured == One,
            $"Integer comparison failed:
                {asigned} > {bsigned} = {res} != {resMeasured} [n={n}]");
        ResetAll(aqs + bqs + [result]);
        for ctrlState in 0..2^numCtrl - 1 {
            ApplyXorInPlace(ctrlState, LittleEndian(ctrlqs));
            ApplyXorInPlace(a, LittleEndian(aqs));
            ApplyXorInPlace(b, LittleEndian(bqs));
            if (signed) {
                (Controlled CompareGTSI) (ctrlqs,
                    (SignedLittleEndian(LittleEndian(aqs)),
                        SignedLittleEndian(LittleEndian(bqs)),
                        result));
            }
            else {
                (Controlled CompareGTI) (ctrlqs,
                    (LittleEndian(aqs),
                        LittleEndian(bqs),
                        result));
            }
            set res = asigned > bsigned;
            if (ctrlState < 2^numCtrl-1) {
                set res = false;
            }
            set resMeasured = M(result);
            EqualityFactB(res, resMeasured == One,
                $"Controlled integer comparison failed.");
            ResetAll(aqs + bqs + [result] + ctrlqs);
        }
    }

    internal operation IntegerMultiplicationRun(signed: Bool, a: Int, b: Int,
                                                na: Int, nb : Int, numCtrl: Int) : Unit {
        let nc = na + nb;
        use aqs = Qubit[na];
        use bqs = Qubit[nb];
        use cqs = Qubit[nc];
        use ctrlqs = Qubit[numCtrl];

        let aLE = LittleEndian(aqs);
        let bLE = LittleEndian(bqs);
        let cLE = LittleEndian(cqs);

        ApplyXorInPlace(a, aLE);
        ApplyXorInPlace(b, bLE);
        if signed {
            MultiplySI(SignedLittleEndian(aLE),
                       SignedLittleEndian(bLE),
                       SignedLittleEndian(cLE));
        } else {
            MultiplyI(aLE, bLE, cLE);
        }

        let asigned = signed and a >= 2^(na - 1) ? -2^na + a | a;
        let bsigned = signed and b >= 2^(nb - 1) ? -2^nb + b | b;
        mutable c = asigned * bsigned;
        mutable cMeasured = MeasureInteger(cLE);
        if signed and cMeasured >= 2^(nc-1) {
            set cMeasured = -2^nc + cMeasured;
        }
        EqualityFactI(c, cMeasured,
            $"Multiplication did not yield the correct result:
                {asigned} * {bsigned} = {c} != {cMeasured} [na={na}, nb={nb}]");
        ResetAll(aqs + bqs + cqs);
        for ctrlState in 0..2^numCtrl - 1 {
            ApplyXorInPlace(ctrlState, LittleEndian(ctrlqs));
            ApplyXorInPlace(a, aLE);
            ApplyXorInPlace(b, bLE);
            if signed {
                Controlled MultiplySI(ctrlqs,
                    (SignedLittleEndian(aLE),
                        SignedLittleEndian(bLE),
                        SignedLittleEndian(cLE)));
            } else {
                Controlled MultiplyI(ctrlqs, (aLE, bLE, cLE));
            }
            set c = asigned * bsigned;
            if ctrlState != 2^numCtrl - 1 {
                set c = 0;
            }
            set cMeasured = MeasureInteger(cLE);
            if signed and cMeasured >= 2^(nc - 1) {
                set cMeasured = -2^nc + cMeasured;
            }
            EqualityFactI(c, cMeasured,
                "Controlled multiplication did not yield the correct result.");
            ResetAll(aqs + bqs + cqs + ctrlqs);
        }
    }

    internal operation IntegerSquareRun(signed: Bool, a: Int,
                                        n: Int, numCtrl: Int) : Unit {
        use aqs = Qubit[n];
        use cqs = Qubit[2 * n];
        use ctrlqs = Qubit[numCtrl];
        ApplyXorInPlace(a, LittleEndian(aqs));
        if signed {
            SquareSI(SignedLittleEndian(LittleEndian(aqs)),
                        SignedLittleEndian(LittleEndian(cqs)));
        } else {
            SquareI(LittleEndian(aqs), LittleEndian(cqs));
        }
        mutable signeda = a;
        if signed and a >= 2^(n-1) {
            set signeda = -2^n + a;
        }
        mutable c = signeda * signeda;
        mutable cMeasured = MeasureInteger(LittleEndian(cqs));
        if signed and cMeasured >= 2^(2*n-1) {
            set cMeasured = -2^(2*n) + cMeasured;
        }
        EqualityFactI(c, cMeasured,
            "Square did not yield the correct result.");
        ResetAll(aqs + cqs);
        for ctrlState in 0..2^numCtrl - 1 {
            ApplyXorInPlace(ctrlState, LittleEndian(ctrlqs));
            ApplyXorInPlace(a, LittleEndian(aqs));
            if signed {
                (Controlled SquareSI) (ctrlqs,
                    (SignedLittleEndian(LittleEndian(aqs)),
                        SignedLittleEndian(LittleEndian(cqs))));
            }
            else {
                (Controlled SquareI) (ctrlqs,
                    (LittleEndian(aqs), LittleEndian(cqs)));
            }
            set c = signeda * signeda;
            if (ctrlState != 2^numCtrl-1) {
                set c = 0;
            }
            set cMeasured = MeasureInteger(LittleEndian(cqs));
            if (signed and cMeasured >= 2^(2*n-1)){
                set cMeasured = -2^(2*n) + cMeasured;
            }
            EqualityFactI(c, cMeasured,
                "Controlled square did not yield the correct result.");
            ResetAll(aqs + cqs + ctrlqs);
        }
    }

    internal operation IntegerReciprocalRun(signed: Bool, a: Int,
                                            n: Int, numCtrl: Int) : Unit {
        use aqs = Qubit[n];
        use cqs = Qubit[2 * n];
        use ctrlqs = Qubit[numCtrl];
        ApplyXorInPlace(a, LittleEndian(aqs));
        ComputeReciprocalI(LittleEndian(aqs), LittleEndian(cqs));
        mutable c = a > 0 ? 2^(2*n-1) / a | (2^(2*n)-1);
        mutable cMeasured = MeasureInteger(LittleEndian(cqs));
        mutable aMeasured = MeasureInteger(LittleEndian(aqs));
        EqualityFactI(a, aMeasured,
            "Reciprocal modified the input.");
        EqualityFactI(c, cMeasured,
            $"Reciprocal did not yield the correct result:
                1/{a} = {c} != {cMeasured}.");
        ResetAll(aqs + cqs);
        for ctrlState in 0..2^numCtrl - 1 {
            ApplyXorInPlace(ctrlState, LittleEndian(ctrlqs));
            ApplyXorInPlace(a, LittleEndian(aqs));
            (Controlled ComputeReciprocalI) (ctrlqs,
                (LittleEndian(aqs), LittleEndian(cqs)));
            set c = a > 0 ? 2^(2*n-1) / a | (2^(2*n)-1);
            if (ctrlState != 2^numCtrl-1) {
                set c = 0;
            }
            set cMeasured = MeasureInteger(LittleEndian(cqs));
            EqualityFactI(c, cMeasured,
                $"Controlled reciprocal did not yield the correct result:
                    1/{a} = {c} != {cMeasured}.");
            set aMeasured = MeasureInteger(LittleEndian(aqs));
            EqualityFactI(a, aMeasured,
                "Reciprocal modified the input.");
            ResetAll(aqs + cqs + ctrlqs);
        }
    }

    internal operation IntegerDivisionRun(a: Int, b: Int, n: Int, numCtrl: Int): Unit {
        use aqs = Qubit[n];
        use bqs = Qubit[n];
        use cqs = Qubit[n];
        use ctrlqs = Qubit[numCtrl];
        mutable c = 0;
        if (b > 0) {
            set c = a / b;
        }
        else{
            set c = 2^n - 1;
        }
        let rem = a - c * b;
        ApplyXorInPlace(a, LittleEndian(aqs));
        ApplyXorInPlace(b, LittleEndian(bqs));
        DivideI (LittleEndian(aqs), LittleEndian(bqs),
                            LittleEndian(cqs));
        mutable resMeasured = MeasureInteger(LittleEndian(cqs));
        mutable remMeasured = MeasureInteger(LittleEndian(aqs));
        EqualityFactI(c, resMeasured,
            $"Controlled division did not yield the correct result:
                {a} / {b} = {c} != {resMeasured}");
        EqualityFactI(remMeasured, rem,
            "Controlled division did not yield the correct remainder:
                rem = {rem} != {remMeasured} for {a} / {b}");
        ResetAll(aqs + bqs + cqs + ctrlqs);
        for ctrlState in 0..2^numCtrl - 1 {
            ApplyXorInPlace(ctrlState, LittleEndian(ctrlqs));
            ApplyXorInPlace(a, LittleEndian(aqs));
            ApplyXorInPlace(b, LittleEndian(bqs));
            (Controlled DivideI) (ctrlqs,
                (LittleEndian(aqs), LittleEndian(bqs), LittleEndian(cqs)));
            set resMeasured = MeasureInteger(LittleEndian(cqs));
            set remMeasured = MeasureInteger(LittleEndian(aqs));
            if ctrlState == 2^numCtrl - 1 {
                EqualityFactI(c, resMeasured,
                    $"Controlled division did not yield the correct result:
                        {a} / {b} = {c} != {resMeasured}");
                EqualityFactI(remMeasured, rem,
                    "Controlled division did not yield the correct remainder:
                        rem = {rem} != {remMeasured} for {a} / {b}");
            } else {
                EqualityFactI(0, resMeasured,
                    "Controlled division was not trivial.");
                EqualityFactI(a, remMeasured,
                    "Controlled division was not trivial.");
            }
            ResetAll(aqs + bqs + cqs + ctrlqs);
        }
    }

    internal operation ExhaustiveTestHelper1Arg(TestFunction: (Int, Int, Int) => Unit) : Unit {
        for numCtrlQubits in 0..2 {
            for numQubits in 1..5 {
                for a in 0..2^numQubits - 1 {
                    TestFunction(a, numQubits, numCtrlQubits);
                }
            }
        }
    }

    // Tests an operation that expects two input arguments with the same number of bits
    internal operation ExhaustiveTestHelper2RegularArgs(TestFunction: (Int, Int, Int, Int) => Unit) : Unit {
        for numCtrlQubits in 0..2 {
            for numQubits in 1..5 {
                for a in 0..2^numQubits - 1 {
                    for b in 0..2^numQubits - 1 {
                        TestFunction(a, b, numQubits, numCtrlQubits);
                    }
                }
            }
        }
    }

    // Tests an operation that expects two input arguments with a different number of bits
    internal operation ExhaustiveTestHelper2NonRegularArgs(TestFunction: (Int, Int, Int, Int, Int) => Unit) : Unit {
        for numCtrlQubits in 0..2 {
            for numQubitsA in 1..4 {
                for numQubitsB in 1..4 {
                    for a in 0..2^numQubitsA - 1 {
                        for b in 0..2^numQubitsB - 1 {
                            TestFunction(a, b, numQubitsA, numQubitsB, numCtrlQubits);
                        }
                    }
                }
            }
        }
    }
}
