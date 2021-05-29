// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Arithmetic {
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Arrays;

    /// # Summary
    /// Automatically chooses between addition with
    /// carry and without, depending on the register size of `ys`.
    ///
    /// # Input
    /// ## xs
    /// $n$-bit addend.
    /// ## ys
    /// Addend with at least $n$ qubits. Will hold the result.
    operation AddI (xs: LittleEndian, ys: LittleEndian) : Unit is Adj + Ctl {
        if Length(xs!) == Length(ys!) {
            RippleCarryAdderNoCarryTTK(xs, ys);
        }
        elif Length(ys!) > Length(xs!) {
            use qs = Qubit[Length(ys!) - Length(xs!) - 1];
            RippleCarryAdderTTK(LittleEndian(xs! + qs),
                                LittleEndian(Most(ys!)), Tail(ys!));
        }
        else {
            fail "xs must not contain more qubits than ys!";
        }
    }

    /// # Summary
    /// Wrapper for integer comparison: `result = x > y`.
    ///
    /// # Input
    /// ## xs
    /// First $n$-bit number
    /// ## ys
    /// Second $n$-bit number
    /// ## result
    /// Will be flipped if $x > y$
    operation CompareGTI (xs: LittleEndian, ys: LittleEndian,
                            result: Qubit) : Unit is Adj + Ctl {
        GreaterThan(xs, ys, result);
    }

    /// # Summary
    /// Wrapper for signed integer comparison: `result = xs > ys`.
    ///
    /// # Input
    /// ## xs
    /// First $n$-bit number
    /// ## ys
    /// Second $n$-bit number
    /// ## result
    /// Will be flipped if $xs > ys$
    operation CompareGTSI (xs: SignedLittleEndian,
                           ys: SignedLittleEndian,
                           result: Qubit) : Unit is Adj + Ctl {
        use tmp = Qubit();
        CNOT(Tail(xs!!), tmp);
        CNOT(Tail(ys!!), tmp);
        X(tmp);
        Controlled CompareGTI([tmp], (xs!, ys!, result));
        X(tmp);
        CCNOT(tmp, Tail(ys!!), result);
        CNOT(Tail(xs!!), tmp);
        CNOT(Tail(ys!!), tmp);
    }
}
