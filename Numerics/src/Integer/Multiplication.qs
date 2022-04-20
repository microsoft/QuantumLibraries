// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Arithmetic {
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Diagnostics;

    /// # Summary
    /// Multiply integer `xs` by integer `ys` and store the result in `result`,
    /// which must be zero initially.
    ///
    /// # Input
    /// ## xs
    /// 𝑛₁-bit multiplicand
    /// ## ys
    /// 𝑛₂-bit multiplier
    /// ## result
    /// (𝑛₁+𝑛₂)-bit result, must be in state |0⟩ initially.
    ///
    /// # Remarks
    /// Uses a standard shift-and-add approach to implement the multiplication.
    /// The controlled version was improved by copying out 𝑥ᵢ to an ancilla
    /// qubit conditioned on the control qubits, and then controlling the
    /// addition on the ancilla qubit.
    operation MultiplyI(xs: LittleEndian, ys: LittleEndian, result: LittleEndian) : Unit is Adj + Ctl {
        body (...) {
            let na = Length(xs!);
            let nb = Length(ys!);

            EqualityFactI(na + nb, Length(result!), "Integer multiplication requires a register as long as both input registers added");
            AssertAllZero(result!);

            for (idx, actl) in Enumerated(xs!) {
                Controlled AddI([actl], (ys, LittleEndian(result![idx..idx + nb])));
            }
        }
        controlled (controls, ...) {
            let na = Length(xs!);
            let nb = Length(ys!);

            EqualityFactI(na + nb, Length(result!), "Integer multiplication requires a register as long as both input registers added");
            AssertAllZero(result!);

            // Perform various optimizations based on number of controls
            let numControls = Length(controls);
            if numControls == 0 {
                MultiplyI(xs, ys, result);
            } elif numControls == 1 {
                use aux = Qubit();
                for (idx, actl) in Enumerated(xs!) {
                    within {
                        ApplyAnd(controls[0], actl, aux);
                    } apply {
                        Controlled AddI([aux], (ys, LittleEndian(result![idx..idx + nb])));
                    }
                }
            } else {
                use helper = Qubit[numControls];
                within {
                    AndLadder(CCNOTop(ApplyAnd), controls, Most(helper));
                } apply {
                    for (idx, actl) in Enumerated(xs!) {
                        within {
                            ApplyAnd(Tail(Most(helper)), actl, Tail(helper));
                        } apply {
                            Controlled AddI([Tail(helper)], (ys, LittleEndian(result![idx..idx + nb])));
                        }
                    }
                }
            }
        }
    }

    /// # Summary
    /// Computes the square of the integer `xs` into `result`,
    /// which must be zero initially.
    ///
    /// # Input
    /// ## xs
    /// 𝑛-bit number to square
    /// ## result
    /// 2𝑛-bit result, must be in state |0⟩ initially.
    ///
    /// # Remarks
    /// Uses a standard shift-and-add approach to compute the square. Saves
    /// 𝑛-1 qubits compared to the straight-forward solution which first
    /// copies out `xs` before applying a regular multiplier and then undoing
    /// the copy operation.
    operation SquareI(xs: LittleEndian, result: LittleEndian) : Unit {
        body (...) {
            Controlled SquareI([], (xs, result));
        }
        controlled (controls, ...) {
            let n = Length(xs!);

            EqualityFactI(2 * n, Length(result!), "Integer multiplication requires a 2n-bit result registers.");
            AssertAllZero(result!);

            let numControls = Length(controls);
            if numControls == 0 {
                use aux = Qubit();
                for (idx, ctl) in Enumerated(xs!) {
                    within {
                        CNOT(ctl, aux);
                    } apply {
                        Controlled AddI([aux], (xs, LittleEndian(result![idx..idx + n])));
                    }
                }
            } elif numControls == 1 {
                use aux = Qubit();
                for (idx, ctl) in Enumerated(xs!) {
                    within {
                        ApplyAnd(controls[0], ctl, aux);
                    } apply {
                        Controlled AddI([aux], (xs, LittleEndian(result![idx..idx + n])));
                    }
                }
            } else {
                use helper = Qubit[numControls];
                within {
                    AndLadder(CCNOTop(ApplyAnd), controls, Most(helper));
                } apply {
                    for (idx, ctl) in Enumerated(xs!) {
                        within {
                            ApplyAnd(Tail(Most(helper)), ctl, Tail(helper));
                        } apply {
                            Controlled AddI([Tail(helper)], (xs, LittleEndian(result![idx..idx + n])));
                        }
                    }
                }
            }
        }
        adjoint auto;
        adjoint controlled auto;
    }

    /// # Summary
    /// Multiply signed integer `xs` by signed integer `ys` and store
    /// the result in `result`, which must be zero initially.
    ///
    /// # Input
    /// ## xs
    /// 𝑛₁-bit multiplicand
    /// ## ys
    /// 𝑛₂-bit multiplier
    /// ## result
    /// (𝑛₁+𝑛₂)-bit result, must be in state |0⟩
    /// initially.
    operation MultiplySI(xs: SignedLittleEndian, ys: SignedLittleEndian, result: SignedLittleEndian): Unit {
        body (...) {
            Controlled MultiplySI([], (xs, ys, result));
        }
        controlled (controls, ...) {
            use signx = Qubit();
            use signy = Qubit();

            within {
                CNOT(Tail(xs!!), signx);
                CNOT(Tail(ys!!), signy);
                Controlled Invert2sSI([signx], xs);
                Controlled Invert2sSI([signy], ys);
            } apply {
                Controlled MultiplyI(controls, (xs!, ys!, result!));
                within {
                    CNOT(signx, signy);
                } apply {
                    // No controls required since `result` will still be zero
                    // if we did not perform the multiplication above.
                    Controlled Invert2sSI([signy], result);
                }
            }
        }
        adjoint auto;
        adjoint controlled auto;
    }

    /// # Summary
    /// Square signed integer `xs` and store
    /// the result in `result`, which must be zero initially.
    ///
    /// # Input
    /// ## xs
    /// 𝑛-bit integer to square
    /// ## result
    /// 2𝑛-bit result, must be in state |0⟩
    /// initially.
    ///
    /// # Remarks
    /// The implementation relies on `SquareI`.
    operation SquareSI (xs: SignedLittleEndian, result: SignedLittleEndian): Unit is Adj + Ctl {
        body (...) {
            Controlled SquareSI([], (xs, result));
        }
        controlled (controls, ...) {
            let n = Length(xs!!);
            use signx = Qubit();
            use signy = Qubit();

            within {
                CNOT(Tail(xs!!), signx);
                Controlled Invert2sSI([signx], xs);
            } apply {
                Controlled SquareI(controls, (xs!, result!));
            }
        }
    }
}
