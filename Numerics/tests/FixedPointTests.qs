// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Numerics.ToffoliTests {
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Convert;
    open Microsoft.Quantum.Math;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Arithmetic;
    open Microsoft.Quantum.Diagnostics;

    operation PrepareFxPTest() : Unit {
        for a in [1.2, 3.9, 3.14159, -0.6, -4.5, -3.1931, 0.0] {
            use xs = Qubit[10];
            let fp = FixedPoint(4, xs);
            PrepareFxP(a, fp);
            let measured = MeasureFxP(fp);
            Fact(AbsD(measured - a) <= 1./IntAsDouble(2^7),
                $"FixedPoint initialized to {a} but measured {measured}.");
            ResetAll(xs);
        }
    }

    internal operation PrepareAsSignedAndMeasure(value : Int, fxp : FixedPoint) : Double {
        ApplyXorInPlace(value, LittleEndian(Snd(fxp!)));
        return MeasureFxP(fxp);
    }

    operation MeasureFxPTest() : Unit {
        use qs = Qubit[4];
        let qsFxP = FixedPoint(2, qs);

        NearEqualityFactD(PrepareAsSignedAndMeasure(0b0000, qsFxP), 0.0);
        NearEqualityFactD(PrepareAsSignedAndMeasure(0b0001, qsFxP), 0.25);
        NearEqualityFactD(PrepareAsSignedAndMeasure(0b0010, qsFxP), 0.5);
        NearEqualityFactD(PrepareAsSignedAndMeasure(0b0011, qsFxP), 0.75);
        NearEqualityFactD(PrepareAsSignedAndMeasure(0b0100, qsFxP), 1.0);
        NearEqualityFactD(PrepareAsSignedAndMeasure(0b0101, qsFxP), 1.25);
        NearEqualityFactD(PrepareAsSignedAndMeasure(0b0110, qsFxP), 1.5);
        NearEqualityFactD(PrepareAsSignedAndMeasure(0b0111, qsFxP), 1.75);
        NearEqualityFactD(PrepareAsSignedAndMeasure(0b1000, qsFxP), -2.0);
        NearEqualityFactD(PrepareAsSignedAndMeasure(0b1001, qsFxP), -1.75);
        NearEqualityFactD(PrepareAsSignedAndMeasure(0b1010, qsFxP), -1.5);
        NearEqualityFactD(PrepareAsSignedAndMeasure(0b1011, qsFxP), -1.25);
        NearEqualityFactD(PrepareAsSignedAndMeasure(0b1100, qsFxP), -1.00);
        NearEqualityFactD(PrepareAsSignedAndMeasure(0b1101, qsFxP), -0.75);
        NearEqualityFactD(PrepareAsSignedAndMeasure(0b1110, qsFxP), -0.5);
        NearEqualityFactD(PrepareAsSignedAndMeasure(0b1111, qsFxP), -0.25);
    }

    operation CompareGreaterThanFxPTest() : Unit {
        for a in [1.2, 3.9, 3.14159, -0.6, -4.5, -3.1931, 0.0] {
            for b in [1.1, 3.95, 3.14259, -0.4, -4.6, -3.931, 0.1] {
                use xs = Qubit[10];
                use ys = Qubit[10];
                use res = Qubit();

                let fp1 = FixedPoint(4, xs);
                let fp2 = FixedPoint(4, ys);
                PrepareFxP(a, fp1);
                PrepareFxP(b, fp2);
                CompareGreaterThanFxP(fp1, fp2, res);
                let measured = M(res);
                EqualityFactB(a > b, measured == One,
                    $"FixedPoint comparison: {a} > {b} != {measured}.");
                ResetAll(xs + ys + [res]);
            }
        }
    }

    operation AddConstantFxPTest() : Unit {
        for a in [1.2, 3.9, 3.14159, -0.6, -4.5, -3.1931, 0.0] {
            for b in [1.2, 3.9, 3.14159, -0.6, -4.5, -3.1931, 0.0] {
                use xs = Qubit[11];
                let fp = FixedPoint(5, xs);
                PrepareFxP(a, fp);
                AddConstantFxP(b, fp);
                let measured = MeasureFxP(fp);
                EqualityWithinToleranceFact(measured, (a+b), 1. / IntAsDouble(2^6));
            }
        }
    }

    operation AddFxPTest() : Unit {
        for a in [1.2, 3.9, 3.14159, -0.6, -4.5, -3.1931, 0.0] {
            for b in [1.2, 3.9, 3.14159, -0.6, -4.5, -3.1931, 0.0] {
                use xs = Qubit[11];
                use ys = Qubit[11];
                let fp1 = FixedPoint(5, xs);
                let fp2 = FixedPoint(5, ys);
                PrepareFxP(a, fp1);
                PrepareFxP(b, fp2);
                AddFxP(fp1, fp2);
                let measured = MeasureFxP(fp2);
                EqualityWithinToleranceFact(measured, (a+b), 1. / IntAsDouble(2^6));
                ResetAll(xs + ys);
            }
        }
    }

    operation MultiplyFxPTest() : Unit {
        for pos in 5..8 {
            for a in [1.2, 3.9, 3.14159, -0.6, -3.5, -3.1931, 0.0] {
                for b in [1.2, 3.9, 3.14159, -0.6, -3.5, -3.1931, 0.0] {
                    use xs = Qubit[13];
                    use ys = Qubit[13];
                    use zs = Qubit[13];

                    let fp1 = FixedPoint(pos, xs);
                    let fp2 = FixedPoint(pos, ys);
                    let fp3 = FixedPoint(pos, zs);
                    PrepareFxP(a, fp1);
                    PrepareFxP(b, fp2);
                    MultiplyFxP(fp1, fp2, fp3);
                    let measured = MeasureFxP(fp3);
                    let eps = 1./IntAsDouble(2^(13-pos));
                    let epsTotal = AbsD(a) * eps + AbsD(b) * eps + eps * eps;
                    EqualityWithinToleranceFact(measured, a * b, epsTotal);
                    ResetAll(xs + ys + zs);
                }
            }
        }
    }

    operation SquareFxPTest() : Unit {
        for pos in 5..8 {
            for a in [1.2, 3.9, 3.14159, -0.6, -3.5, -3.1931, 0.0] {
                use xs = Qubit[13];
                use ys = Qubit[13];
                let fp1 = FixedPoint(pos, xs);
                let fp2 = FixedPoint(pos, ys);
                PrepareFxP(a, fp1);
                SquareFxP(fp1, fp2);
                let measured = MeasureFxP(fp2);
                let eps = 1./IntAsDouble(2^(13-pos));
                let epsTotal = 2. * AbsD(a) * eps + eps * eps;
                EqualityWithinToleranceFact(measured, a * a, epsTotal);
                ResetAll(xs + ys);
            }
        }
    }

    function _computeReciprocal(a : Double, n : Int, pos : Int, pos2 : Int) : Double {
        let p = pos;
        let intA = a >= 0. ? Floor(AbsD(a) * IntAsDouble(2^(n-p)) + 0.5)
                           | Ceiling(AbsD(a) * IntAsDouble(2^(n-p)) - 0.5);
        let intDiv = 2^(2*n-1) / intA;
        let aReciprUnsigned = IntAsDouble((intDiv >>> (p+pos2-1)) &&& (2^n-1)) / IntAsDouble(2^(n-pos2));
        return (a >= 0. ? 1. | -1.) * aReciprUnsigned;
    }

    operation ComputeReciprocalFxPTest() : Unit {
        for pos in 5..8 {
            for pos2 in pos - 1..pos + 3 {
                for a in [1.2, 3.9, -0.314159, -0.6, -3.5, -3.1931, 0.127] {
                    let n = 20;
                    use xs = Qubit[n];
                    use ys = Qubit[n];
                    let fp1 = FixedPoint(pos, xs);
                    let fp2 = FixedPoint(pos, ys);
                    PrepareFxP(a, fp1);
                    ComputeReciprocalFxP(fp1, fp2);
                    let measured = MeasureFxP(fp2);
                    let eps = 1./IntAsDouble(2^(n-pos));
                    let eps2 = 1./IntAsDouble(2^(n-pos2));
                    let aEpsLarger = a + (a>=0. ? eps | -eps);
                    let aEpsSmaller = a - (a>=0. ? eps | -eps);
                    let res1 = _computeReciprocal(a+eps,n,pos,pos2);
                    let res2 = _computeReciprocal(a-eps,n,pos,pos2);
                    let minRes = MinD(res1, res2) - eps2;
                    let maxRes = MaxD(res1, res2) + eps2;
                    let isWithinTol = minRes <= measured and
                                        maxRes >= measured;
                    Fact(isWithinTol,
                        $"FixedPoint reciprocal 1/{a}: {measured} is not within [{minRes},{maxRes}].");
                    ResetAll(xs + ys);
                }
            }
        }
    }

    operation SquareFxPCtrlTest() : Unit {
        for ctrl in 0..3 {
            for pos in 5..8 {
                for a in [1.2, 3.9, 3.14159, -0.6, -3.5, -3.1931, 0.0] {
                    use xs = Qubit[13];
                    use ys = Qubit[13];
                    use cs = Qubit[2];

                    ApplyXorInPlace(ctrl, LittleEndian(cs));
                    let fp1 = FixedPoint(pos, xs);
                    let fp2 = FixedPoint(pos, ys);
                    PrepareFxP(a, fp1);
                    (Controlled SquareFxP)(cs, (fp1, fp2));
                    let measured = MeasureFxP(fp2);
                    let eps = 1./IntAsDouble(2^(13-pos));
                    let epsTotal = 2. * AbsD(a) * eps + eps * eps;
                    if (ctrl == 3) {
                        EqualityWithinToleranceFact(measured, a * a, epsTotal);
                    }
                    else {
                        let measuredI = MeasureInteger(LittleEndian(ys));
                        EqualityFactI(measuredI, 0,
                            "Controlled FixedPoint square changed the result register!");
                    }
                    ResetAll(xs + ys + cs);
                }
            }
        }
    }

    operation MultiplyFxPCtrlTest() : Unit {
        for ctrl in 0..3 {
            for pos in 5..8 {
                for a in [1.2, 3.9, 3.14159, -0.6, -3.5, -3.1931, 0.0] {
                    for b in [1.2, 3.9, 3.14159, -0.6, -3.5, -3.1931, 0.0] {
                        use xs = Qubit[13];
                        use ys = Qubit[13];
                        use zs = Qubit[13];
                        use cs = Qubit[2];
                        ApplyXorInPlace(ctrl, LittleEndian(cs));
                        let fp1 = FixedPoint(pos, xs);
                        let fp2 = FixedPoint(pos, ys);
                        let fp3 = FixedPoint(pos, zs);
                        PrepareFxP(a, fp1);
                        PrepareFxP(b, fp2);
                        (Controlled MultiplyFxP)(cs, (fp1, fp2, fp3));
                        let measured = MeasureFxP(fp3);
                        let eps = 1./IntAsDouble(2^(13-pos));
                        let epsTotal = AbsD(a) * eps + AbsD(b) * eps + eps * eps;
                        if ctrl == 3 {
                            EqualityWithinToleranceFact(measured, a * b, epsTotal);
                        }
                        else {
                            let measuredI = MeasureInteger(LittleEndian(zs));
                            EqualityFactI(measuredI, 0,
                                "Controlled FixedPoint multiplication changed the output register!");
                        }
                        ResetAll(xs + ys + zs + cs);
                    }
                }
            }
        }
    }

    operation EvaluatePolynomialFxPTest() : Unit {
        for pos in 4..5 {
            for coeffs in [[1.3, -2.4, 1.9],
                            [-0.3, -0.2],
                            [0.1, 1.1, -0.1, 0.2],
                            [0.1, -0.1, 0.2, 0.2, -0.1, 0.3],
                            [0.2]] {
                for a in [0.0, 0.1, -0.1, 0.2, -0.2, 1.3, -1.3] {
                    let n = 20;
                    use xs = Qubit[n];
                    use ys = Qubit[n];
                    let fp1 = FixedPoint(pos, xs);
                    let fp2 = FixedPoint(pos, ys);
                    PrepareFxP(a, fp1);
                    EvaluatePolynomialFxP(coeffs, fp1, fp2);
                    let measured = MeasureFxP(fp2);
                    let eps = 1./IntAsDouble(2^(n-pos));
                    mutable epsTotal = 0.;
                    mutable errX = eps;
                    mutable result = Tail(coeffs);
                    set epsTotal = epsTotal + eps;
                    for coeff in coeffs[Length(coeffs)-2..(-1)..0] {
                        set epsTotal = epsTotal + AbsD(result) * eps
                                        + AbsD(a) * epsTotal + eps * epsTotal;
                        set result = result * a + coeff;
                        set epsTotal += eps;
                    }
                    EqualityWithinToleranceFact(measured, result, epsTotal);
                    ResetAll(xs + ys);
                }
            }
        }
    }

    operation EvaluatePolynomialFxPCtrlTest() : Unit {
        for ctrl in 0..3 {
            for pos in 4..5 {
                for coeffs in [[1.3, -2.4, 1.9],
                                [-0.3, -0.2],
                                [0.2]] {
                    for a in [0.0, 0.1, -0.2, 1.3, -1.3] {
                        let n = 20;
                        use xs = Qubit[n];
                        use ys = Qubit[n];
                        use ctrls = Qubit[2];
                        let fp1 = FixedPoint(pos, xs);
                        let fp2 = FixedPoint(pos, ys);
                        ApplyXorInPlace(ctrl, LittleEndian(ctrls));
                        PrepareFxP(a, fp1);
                        (Controlled EvaluatePolynomialFxP)(ctrls,
                            (coeffs, fp1, fp2));
                        let measured = MeasureFxP(fp2);
                        let eps = 1./IntAsDouble(2^(n-pos));
                        mutable epsTotal = 0.;
                        mutable errX = eps;
                        mutable result = Tail(coeffs);
                        set epsTotal = epsTotal + eps;
                        for coeff in coeffs[Length(coeffs) - 2..(-1)..0] {
                            set epsTotal = epsTotal + AbsD(result) * eps
                                            + AbsD(a) * epsTotal + eps * epsTotal;
                            set result = result * a + coeff;
                            set epsTotal += eps;
                        }
                        if ctrl == 3 {
                            EqualityWithinToleranceFact(measured, result, epsTotal);
                        }
                        else{
                            let measuredI = MeasureInteger(LittleEndian(ys));
                            EqualityFactI(measuredI,
                                0,
                                $"Controlled FixedPoint polynomial evaluation changed the output register!");
                        }
                        ResetAll(xs + ys + ctrls);
                    }
                }
            }
        }
    }

    operation EvaluateOddPolynomialFxPTest() : Unit {
        for pos in 4..5 {
            for coeffs in [[1.3, -2.4, 1.9],
                            [-0.3],
                            [0.1, -0.1, 0.2, 0.2]] {
                for a in [0.0, 0.1, -0.1, 0.2, -0.2, 1.3, -1.3] {
                    let n = 20;
                    use xs = Qubit[n];
                    use ys = Qubit[n];
                    let fp1 = FixedPoint(pos, xs);
                    let fp2 = FixedPoint(pos, ys);
                    PrepareFxP(a, fp1);
                    EvaluateOddPolynomialFxP(coeffs, fp1, fp2);
                    let measured = MeasureFxP(fp2);
                    let eps = 1./IntAsDouble(2^(n-pos));
                    mutable epsTotal = 0.;
                    mutable errX = eps;
                    mutable result = Tail(coeffs);
                    set epsTotal = epsTotal + eps;
                    let aSquare = a * a;
                    for coeff in coeffs[Length(coeffs)-2..(-1)..0] {
                        set epsTotal = epsTotal + AbsD(result) * eps
                                        + AbsD(aSquare) * epsTotal + eps * epsTotal;
                        set result = result * aSquare + coeff;
                        set epsTotal = epsTotal + eps;
                    }
                    set epsTotal = epsTotal + AbsD(result) * eps + AbsD(a) * epsTotal
                                    + eps * epsTotal;
                    set result = result * a;

                    EqualityWithinToleranceFact(measured, result, epsTotal);
                    ResetAll(xs + ys);
                }
            }
        }
    }

    operation EvaluateOddPolynomialFxPCtrlTest() : Unit {
        for ctrl in 0..3 {
            for pos in 4..5 {
                for coeffs in [[1.3, -2.4, 1.9],
                                [-0.3, -0.2],
                                [0.1, 1.1, -0.1, 0.2],
                                [0.1, -0.1, 0.2, 0.2, -0.1, 0.3],
                                [0.2]] {
                    for a in [0.0, 0.1, -0.1, 0.2, -0.2, 1.3, -1.3] {
                        let n = 20;
                        use xs = Qubit[n];
                        use ys = Qubit[n];
                        use ctrls = Qubit[2];
                        let fp1 = FixedPoint(pos, xs);
                        let fp2 = FixedPoint(pos, ys);
                        PrepareFxP(a, fp1);
                        EvaluateOddPolynomialFxP(coeffs, fp1, fp2);
                        let measured = MeasureFxP(fp2);
                        let eps = 1./IntAsDouble(2^(n-pos));
                        mutable epsTotal = 0.;
                        mutable errX = eps;
                        mutable result = Tail(coeffs);
                        set epsTotal = epsTotal + eps;
                        let aSquare = a * a;
                        for coeff in coeffs[Length(coeffs)-2..(-1)..0] {
                            set epsTotal = epsTotal + AbsD(result) * eps
                                            + AbsD(aSquare) * epsTotal + eps * epsTotal;
                            set result = result * aSquare + coeff;
                            set epsTotal = epsTotal + eps;
                        }
                        set epsTotal = epsTotal + AbsD(result) * eps + AbsD(a) * epsTotal
                                        + eps * epsTotal;
                        set result = result * a;
                        if (ctrl == 3) {
                            EqualityWithinToleranceFact(measured, result, epsTotal);
                        }
                        else{
                            let measuredI = MeasureInteger(LittleEndian(ys));
                            EqualityFactI(measuredI,
                                0,
                                $"Controlled FixedPoint polynomial evaluation changed the output register!");
                        }
                        ResetAll(xs + ys + ctrls);
                    }
                }
            }
        }
    }
}
