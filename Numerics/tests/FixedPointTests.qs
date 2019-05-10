// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Numerics.ToffoliTests {
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Convert;
    open Microsoft.Quantum.Math;
    open Microsoft.Quantum.Primitive;
    open Microsoft.Quantum.Arithmetic;
	open Microsoft.Quantum.Diagnostics;

    operation FixedPointInitTest() : Unit {
        for (a in [1.2, 3.9, 3.14159, -0.6, -4.5, -3.1931, 0.0]){
            using (xs = Qubit[10]) {
                let fp = FixedPoint(4, xs);
                FixedPointInit(a, fp);
                let measured = MeasureFixedPoint(fp);
                EqualityFactB(AbsD(measured - a) <= 1./IntAsDouble(2^7), true,
                    $"FixedPoint initialized to {a} but measured {measured}.");
                ResetAll(xs);
            }
        }
    }

    operation FixedPointGreaterThanTest() : Unit {
        for (a in [1.2, 3.9, 3.14159, -0.6, -4.5, -3.1931, 0.0]){
            for (b in [1.1, 3.95, 3.14259, -0.4, -4.6, -3.931, 0.1]) {
                using ((xs, ys, res) = (Qubit[10], Qubit[10], Qubit())) {
                    let fp1 = FixedPoint(4, xs);
                    let fp2 = FixedPoint(4, ys);
                    FixedPointInit(a, fp1);
                    FixedPointInit(b, fp2);
                    FixedPointGreaterThan(fp1, fp2, res);
                    let measured = M(res);
                    EqualityFactB(a > b, measured == One,
                        $"FixedPoint comparison: {a} > {b} != {measured}.");
                    ResetAll(xs + ys + [res]);
                }
            }
        }
    }

    operation FixedPointAdditionConstantTest() : Unit {
        for (a in [1.2, 3.9, 3.14159, -0.6, -4.5, -3.1931, 0.0]){
            for (b in [1.2, 3.9, 3.14159, -0.6, -4.5, -3.1931, 0.0]){
                using (xs = Qubit[11]) {
                    let fp = FixedPoint(5, xs);
                    FixedPointInit(a, fp);
                    FixedPointAdditionConstant(b, fp);
                    let measured = MeasureFixedPoint(fp);
                    EqualityFactB(AbsD(measured - (a+b)) <= 1./IntAsDouble(2^6),
                        true,
                        $"FixedPoint addition {a}+{b} != {measured}.");
                }
            }
        }
    }

    operation FixedPointAdditionTest() : Unit {
        for (a in [1.2, 3.9, 3.14159, -0.6, -4.5, -3.1931, 0.0]){
            for (b in [1.2, 3.9, 3.14159, -0.6, -4.5, -3.1931, 0.0]){
                using ((xs, ys) = (Qubit[11], Qubit[11])) {
                    let fp1 = FixedPoint(5, xs);
                    let fp2 = FixedPoint(5, ys);
                    FixedPointInit(a, fp1);
                    FixedPointInit(b, fp2);
                    FixedPointAddition(fp1, fp2);
                    let measured = MeasureFixedPoint(fp2);
                    EqualityFactB(AbsD(measured - (a+b)) <= 1./IntAsDouble(2^6),
                        true,
                        $"FixedPoint addition {a}+{b} != {measured}.");
                    ResetAll(xs + ys);
                }
            }
        }
    }

    operation FixedPointMultiplicationTest() : Unit {
        for (pos in 5..8) {
            for (a in [1.2, 3.9, 3.14159, -0.6, -3.5, -3.1931, 0.0]){
                for (b in [1.2, 3.9, 3.14159, -0.6, -3.5, -3.1931, 0.0]){
                    using ((xs, ys, zs) = (Qubit[13], Qubit[13], Qubit[13])) {
                        let fp1 = FixedPoint(pos, xs);
                        let fp2 = FixedPoint(pos, ys);
                        let fp3 = FixedPoint(pos, zs);
                        FixedPointInit(a, fp1);
                        FixedPointInit(b, fp2);
                        FixedPointMultiplication(fp1, fp2, fp3);
                        let measured = MeasureFixedPoint(fp3);
                        let eps = 1./IntAsDouble(2^(13-pos));
                        let epsTotal = AbsD(a) * eps + AbsD(b) * eps + eps * eps;
                        EqualityFactB(AbsD(measured - (a*b)) <= epsTotal,
                            true,
                            $"FixedPoint multiplication {a}*{b} != {measured}.");
                        ResetAll(xs + ys + zs);
                    }
                }
            }
        }
    }

    operation FixedPointSquareTest() : Unit {
        for (pos in 5..8) {
            for (a in [1.2, 3.9, 3.14159, -0.6, -3.5, -3.1931, 0.0]){
                using ((xs, ys) = (Qubit[13], Qubit[13])) {
                    let fp1 = FixedPoint(pos, xs);
                    let fp2 = FixedPoint(pos, ys);
                    FixedPointInit(a, fp1);
                    FixedPointSquare(fp1, fp2);
                    let measured = MeasureFixedPoint(fp2);
                    let eps = 1./IntAsDouble(2^(13-pos));
                    let epsTotal = 2. * AbsD(a) * eps + eps * eps;
                    EqualityFactB(AbsD(measured - (a*a)) <= epsTotal,
                        true,
                        $"FixedPoint square {a}^2 != {measured}.");
                    ResetAll(xs + ys);
                }
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

    operation FixedPointReciprocalTest() : Unit {
        for (pos in 5..8) {
            for (pos2 in pos-1..pos+3) {
                for (a in [1.2, 3.9, -0.314159, -0.6, -3.5, -3.1931, 0.127]){
                    let n = 20;
                    using ((xs, ys) = (Qubit[n], Qubit[n])) {
                        let fp1 = FixedPoint(pos, xs);
                        let fp2 = FixedPoint(pos, ys);
                        FixedPointInit(a, fp1);
                        FixedPointReciprocal(fp1, fp2);
                        let measured = MeasureFixedPoint(fp2);
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
                        EqualityFactB(isWithinTol,
                            true,
                            $"FixedPoint reciprocal 1/{a}: {measured} is not within [{minRes},{maxRes}].");
                        ResetAll(xs + ys);
                    }
                }
            }
        }
    }

    operation FixedPointSquareCtrlTest() : Unit {
        for (ctrl in 0..3) {
            for (pos in 5..8) {
                for (a in [1.2, 3.9, 3.14159, -0.6, -3.5, -3.1931, 0.0]){
                    using ((xs, ys, cs) = (Qubit[13], Qubit[13], Qubit[2])) {
                        ApplyXorInPlace(ctrl, LittleEndian(cs));
                        let fp1 = FixedPoint(pos, xs);
                        let fp2 = FixedPoint(pos, ys);
                        FixedPointInit(a, fp1);
                        (Controlled FixedPointSquare)(cs, (fp1, fp2));
                        let measured = MeasureFixedPoint(fp2);
                        let eps = 1./IntAsDouble(2^(13-pos));
                        let epsTotal = 2. * AbsD(a) * eps + eps * eps;
                        if (ctrl == 3){
                            EqualityFactB(AbsD(measured - (a*a)) <= epsTotal,
                                true,
                                $"FixedPoint square {a}^2 != {measured}.");
                        }
                        else{
                            let measuredI = MeasureInteger(LittleEndian(ys));
                            EqualityFactI(measuredI, 0,
                                "Controlled FixedPoint square changed the result register!");
                        }
                        ResetAll(xs + ys + cs);
                    }
                }
            }
        }
    }

    operation FixedPointMultiplicationCtrlTest() : Unit {
        for (ctrl in 0..3) {
            for (pos in 5..8) {
                for (a in [1.2, 3.9, 3.14159, -0.6, -3.5, -3.1931, 0.0]){
                    for (b in [1.2, 3.9, 3.14159, -0.6, -3.5, -3.1931, 0.0]){
                        using ((xs, ys, zs, cs) = (Qubit[13], Qubit[13], Qubit[13], Qubit[2])) {
                            ApplyXorInPlace(ctrl, LittleEndian(cs));
                            let fp1 = FixedPoint(pos, xs);
                            let fp2 = FixedPoint(pos, ys);
                            let fp3 = FixedPoint(pos, zs);
                            FixedPointInit(a, fp1);
                            FixedPointInit(b, fp2);
                            (Controlled FixedPointMultiplication)(cs, (fp1, fp2, fp3));
                            let measured = MeasureFixedPoint(fp3);
                            let eps = 1./IntAsDouble(2^(13-pos));
                            let epsTotal = AbsD(a) * eps + AbsD(b) * eps + eps * eps;
                            if (ctrl == 3) {
                                EqualityFactB(AbsD(measured - (a*b)) <= epsTotal,
                                    true,
                                    $"FixedPoint multiplication {a}*{b} != {measured}.");
                            }
                            else{
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
    }

    operation FixedPointPolynomialTest() : Unit {
        for (pos in 4..5) {
            for (coeffs in [[1.3, -2.4, 1.9],
                            [-0.3, -0.2],
                            [0.1, 1.1, -0.1, 0.2],
                            [0.1, -0.1, 0.2, 0.2, -0.1, 0.3],
                            [0.2]]){
                for (a in [0.0, 0.1, -0.1, 0.2, -0.2, 1.3, -1.3]){
                    let n = 20;
                    using ((xs, ys) = (Qubit[n], Qubit[n])) {
                        let fp1 = FixedPoint(pos, xs);
                        let fp2 = FixedPoint(pos, ys);
                        FixedPointInit(a, fp1);
                        FixedPointPolynomial(coeffs, fp1, fp2);
                        let measured = MeasureFixedPoint(fp2);
                        let eps = 1./IntAsDouble(2^(n-pos));
                        mutable epsTotal = 0.;
                        mutable errX = eps;
                        mutable result = Tail(coeffs);
                        set epsTotal = epsTotal + eps;
                        for (coeff in coeffs[Length(coeffs)-2..(-1)..0]) {
                            set epsTotal = epsTotal + AbsD(result) * eps
                                           + AbsD(a) * epsTotal + eps * epsTotal;
                            set result = result * a + coeff;
                            set epsTotal = epsTotal + eps;
                        }
                        EqualityFactB(AbsD(measured - result) <= epsTotal,
                            true,
                            $"FixedPoint polynomial evaluation: {measured} != {result}.");
                        ResetAll(xs + ys);
                    }
                }
            }
        }
    }

    operation FixedPointPolynomialCtrlTest() : Unit {
        for (ctrl in 0..3) {
            for (pos in 4..5) {
                for (coeffs in [[1.3, -2.4, 1.9],
                                [-0.3, -0.2],
                                [0.2]]){
                    for (a in [0.0, 0.1, -0.2, 1.3, -1.3]){
                        let n = 20;
                        using ((xs, ys, ctrls) = (Qubit[n], Qubit[n], Qubit[2])) {
                            let fp1 = FixedPoint(pos, xs);
                            let fp2 = FixedPoint(pos, ys);
                            ApplyXorInPlace(ctrl, LittleEndian(ctrls));
                            FixedPointInit(a, fp1);
                            (Controlled FixedPointPolynomial)(ctrls,
                                (coeffs, fp1, fp2));
                            let measured = MeasureFixedPoint(fp2);
                            let eps = 1./IntAsDouble(2^(n-pos));
                            mutable epsTotal = 0.;
                            mutable errX = eps;
                            mutable result = Tail(coeffs);
                            set epsTotal = epsTotal + eps;
                            for (coeff in coeffs[Length(coeffs)-2..(-1)..0]) {
                                set epsTotal = epsTotal + AbsD(result) * eps
                                               + AbsD(a) * epsTotal + eps * epsTotal;
                                set result = result * a + coeff;
                                set epsTotal = epsTotal + eps;
                            }
                            if (ctrl == 3) {
                                EqualityFactB(AbsD(measured - result) <= epsTotal,
                                    true,
                                    $"FixedPoint polynomial evaluation: {measured} != {result}.");
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

    operation FixedPointOddPolynomialTest() : Unit {
        for (pos in 4..5) {
            for (coeffs in [[1.3, -2.4, 1.9],
                            [-0.3],
                            [0.1, -0.1, 0.2, 0.2]]){
                for (a in [0.0, 0.1, -0.1, 0.2, -0.2, 1.3, -1.3]){
                    let n = 20;
                    using ((xs, ys) = (Qubit[n], Qubit[n])) {
                        let fp1 = FixedPoint(pos, xs);
                        let fp2 = FixedPoint(pos, ys);
                        FixedPointInit(a, fp1);
                        FixedPointOddPolynomial(coeffs, fp1, fp2);
                        let measured = MeasureFixedPoint(fp2);
                        let eps = 1./IntAsDouble(2^(n-pos));
                        mutable epsTotal = 0.;
                        mutable errX = eps;
                        mutable result = Tail(coeffs);
                        set epsTotal = epsTotal + eps;
                        let aSquare = a * a;
                        for (coeff in coeffs[Length(coeffs)-2..(-1)..0]) {
                            set epsTotal = epsTotal + AbsD(result) * eps
                                           + AbsD(aSquare) * epsTotal + eps * epsTotal;
                            set result = result * aSquare + coeff;
                            set epsTotal = epsTotal + eps;
                        }
                        set epsTotal = epsTotal + AbsD(result) * eps + AbsD(a) * epsTotal
                                       + eps * epsTotal;
                        set result = result * a;

                        EqualityFactB(AbsD(measured - result) <= epsTotal,
                            true,
                            $"FixedPoint odd polynomial evaluation: {measured} != {result}.");
                        ResetAll(xs + ys);
                    }
                }
            }
        }
    }

    operation FixedPointOddPolynomialCtrlTest() : Unit {
        for (ctrl in 0..3) {
            for (pos in 4..5) {
                for (coeffs in [[1.3, -2.4, 1.9],
                                [-0.3, -0.2],
                                [0.1, 1.1, -0.1, 0.2],
                                [0.1, -0.1, 0.2, 0.2, -0.1, 0.3],
                                [0.2]]){
                    for (a in [0.0, 0.1, -0.1, 0.2, -0.2, 1.3, -1.3]){
                        let n = 20;
                        using ((xs, ys, ctrls) = (Qubit[n], Qubit[n], Qubit[2])) {
                            let fp1 = FixedPoint(pos, xs);
                            let fp2 = FixedPoint(pos, ys);
                            FixedPointInit(a, fp1);
                            FixedPointOddPolynomial(coeffs, fp1, fp2);
                            let measured = MeasureFixedPoint(fp2);
                            let eps = 1./IntAsDouble(2^(n-pos));
                            mutable epsTotal = 0.;
                            mutable errX = eps;
                            mutable result = Tail(coeffs);
                            set epsTotal = epsTotal + eps;
                            let aSquare = a * a;
                            for (coeff in coeffs[Length(coeffs)-2..(-1)..0]) {
                                set epsTotal = epsTotal + AbsD(result) * eps
                                               + AbsD(aSquare) * epsTotal + eps * epsTotal;
                                set result = result * aSquare + coeff;
                                set epsTotal = epsTotal + eps;
                            }
                            set epsTotal = epsTotal + AbsD(result) * eps + AbsD(a) * epsTotal
                                           + eps * epsTotal;
                            set result = result * a;
                            if (ctrl == 3) {
                                EqualityFactB(AbsD(measured - result) <= epsTotal,
                                    true,
                                    $"FixedPoint polynomial evaluation: {measured} != {result}.");
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
}