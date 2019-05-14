// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Arithmetic {
	open Microsoft.Quantum.Intrinsic;
	open Microsoft.Quantum.Arrays;

    /// # Summary
    /// Wrapper for addition: Automatically chooses between addition with
    /// carry and without, depending on the register size of `ys`.
    ///
    /// # Input
    /// ## xs
    /// $n$-bit addend (LittleEndian)
    /// ## ys
    /// Addend with at least $n$ qubits (LittleEndian). Will hold the result.
    operation IntegerAddition (xs: LittleEndian, ys: LittleEndian) : Unit {
        body (...) {
            if (Length(xs!) == Length(ys!)) {
                RippleCarryAdderNoCarryTTK(xs, ys);
            }
            elif (Length(ys!) > Length(xs!)) {
                using (qs = Qubit[Length(ys!) - Length(xs!) - 1]){
                    RippleCarryAdderTTK(LittleEndian(xs! + qs),
                                        LittleEndian(Most(ys!)), Tail(ys!));
                }
            }
            else{
                fail "xs must not contain more qubits than ys!";
            }
        }
        controlled auto;
        adjoint auto;
        adjoint controlled auto;
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
    operation IntegerGreaterThan (xs: LittleEndian, ys: LittleEndian,
                                  result: Qubit) : Unit {
        body (...) {
            GreaterThan(xs, ys, result);
        }
        controlled auto;
        adjoint auto;
        adjoint controlled auto;
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
    operation SignedIntegerGreaterThan (xs: SignedLittleEndian,
                                        ys: SignedLittleEndian,
                                        result: Qubit) : Unit {
        body (...) {
            using (tmp = Qubit()) {
                CNOT(Tail(xs!!), tmp);
                CNOT(Tail(ys!!), tmp);
                X(tmp);
                (Controlled IntegerGreaterThan)([tmp], (xs!, ys!, result));
                X(tmp);
                CCNOT(tmp, Tail(ys!!), result);
                CNOT(Tail(xs!!), tmp);
                CNOT(Tail(ys!!), tmp);
            }
        }
        controlled auto;
        adjoint auto;
        adjoint controlled auto;
    }
}