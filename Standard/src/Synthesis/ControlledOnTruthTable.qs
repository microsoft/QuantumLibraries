// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Synthesis {
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Convert;
    open Microsoft.Quantum.Diagnostics;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Math;

    /// # Summary
    /// Computes Hadamard transform of a Boolean function in {-1,1} encoding
    /// using Yates's method
    ///
    /// # Input
    /// ## func
    /// Truth table in {-1,1} encoding
    ///
    /// # Output
    /// Spectral coefficients of the function
    ///
    /// # Example
    /// ```Q#
    /// FastHadamardTransform([1, 1, 1, -1]); // [2, 2, 2, -2]
    /// ```
    ///
    /// # Reference
    /// Frank Yates: The design and analysis of factorial experiments, in:
    /// Technical Communication No. 35, Imperial Bureau of Soil Science,
    /// London (1937)
    internal function FastHadamardTransform(func : Int[]) : Int[] {
        let bits = BitSizeI(Length(func) - 1);
        mutable res = func;
        for (m in 0..bits - 1) {
            mutable s = 1 <<< m;
            for (i in 0..(2 * s)..Length(func) - 1) {
                mutable k = i + s;
                for (j in i..i + s - 1) {
                    mutable t = res[j];
                    set res w/= j <- res[j] + res[k];
                    set res w/= k <- t - res[k];
                    set k = k + 1;
                }
            }
        }
        return res;
    }

    /// # Summary
    /// Extends a spectrum by inverted coefficients
    ///
    /// # Input
    /// ## spectrum
    /// Spectral coefficients
    ///
    /// # Output
    /// Coefficients followed by inverted copy
    ///
    /// # Example
    /// ```Q#
    /// Extend([2, 2, 2, -2]); // [2, 2, 2, -2, -2, -2, -2, 2]
    /// ```
    internal function Extend(spectrum : Int[]) : Int[] {
        return spectrum + Mapped(NegationI, spectrum);
    }

    /// # Summary
    /// {-1,1} coding of a Boolean truth value
    ///
    /// # Input
    /// ## b
    /// Boolean value
    ///
    /// # Output
    /// 1, if `b` is false, otherwise -1
    internal function RMEncoding(b : Bool) : Int {
        return b ? -1 | 1;
    }

    /// # Summary
    /// Encode truth table in {1,-1} coding
    ///
    /// # Input
    /// ## table
    /// Truth table as array of truth values
    ///
    /// # Output
    /// Truth table as array of {1,-1} integers
    ///
    /// # Example
    /// ```Q#
    /// Encode([false, false, false, true]); // [1, 1, 1, -1]
    /// ```
    internal function Encode(table : Bool[]) : Int[] {
        return Mapped(RMEncoding, table);
    }

    internal function SizeAdjustedTruthTable(table : Bool[], numVars : Int) : Bool[] {
        let numEntries = 2^numVars;
        if (numEntries < Length(table)) {
            return table[...numEntries - 1];
        } elif (numEntries > Length(table)) {
            return Padded(numEntries, false, table);
        } else {
            return table;
        }
    }

    operation ControlledXOnTruthTable (func : BigInt, controlRegister : Qubit[], targetRegister : Qubit) : Unit {
        body (...) {
            let vars = Length(controlRegister);

            let maxValue = 1L <<< 2^vars;
            Fact(func >= 0L and func < maxValue, $"Argument func must be value from 0 to {maxValue}");

            let tt = BigIntAsBoolArray(func);
            let table = Encode(SizeAdjustedTruthTable(BigIntAsBoolArray(func), vars));
            let spectrum = Extend(FastHadamardTransform(table));

            let qubits = controlRegister + [targetRegister];

            HY(targetRegister);

            for (i in 0..vars) {
                let start = 1 <<< i;
                let code = GrayCode(i);
                for (j in 0..Length(code) - 1) {
                    let (offset, ctrl) = code[j];
                    R1Frac(spectrum[start + offset], vars + 1, qubits[i]);
                    if (i != 0) {
                        CNOT(qubits[ctrl], qubits[i]);
                    }
                }
            }

            H(targetRegister);
        }
        adjoint self;
        controlled (controls, ...) {
            using (q = Qubit()) {
                within {
                    ControlledXOnTruthTableWithCleanTarget(func, controlRegister, q);
                } apply {
                    Controlled X(controls + [q], targetRegister);
                }
            }
        }
        controlled adjoint self;
    }

    operation ControlledXOnTruthTableWithCleanTarget (func : BigInt, controlRegister : Qubit[], targetRegister : Qubit) : Unit {
        body (...) {
            let vars = Length(controlRegister);

            let maxValue = PowL(2L, 2^vars);
            Fact(func >= 0L and func < maxValue, $"Argument func must be value from 0 to {maxValue}");
            AssertAllZero([targetRegister]);

            let tt = BigIntAsBoolArray(func);
            let table = Encode(SizeAdjustedTruthTable(BigIntAsBoolArray(func), vars));
            let spectrum = FastHadamardTransform(table);

            HY(targetRegister);

            let code = GrayCode(vars);
            for (j in 0..Length(code) - 1) {
                let (offset, ctrl) = code[j];
                R1Frac(-spectrum[offset], vars + 1, targetRegister);
                CNOT(controlRegister[ctrl], targetRegister);
            }

            H(targetRegister);
        }
        adjoint (...) {
            let vars = Length(controlRegister);

            let maxValue = PowL(2L, 2^vars);
            Fact(func >= 0L and func < maxValue, $"Argument func must be value from 0 to {maxValue}");

            let tt = BigIntAsBoolArray(func);
            let table = Encode(SizeAdjustedTruthTable(BigIntAsBoolArray(func), vars));
            let spectrum = FastHadamardTransform(table);

            H(targetRegister);
            AssertProb([PauliZ], [targetRegister], One, 0.5, "Probability of the measurement must be 0.5", 1e-10);

            if (IsResultOne(M(targetRegister))) {
                for (i in 0..vars - 1) {
                    let start = 1 <<< i;
                    let code = GrayCode(i);
                    for (j in 0..Length(code) - 1) {
                        let (offset, ctrl) = code[j];
                        R1Frac(spectrum[start + offset], vars, controlRegister[i]);
                        if (i != 0) {
                            CNOT(controlRegister[ctrl], controlRegister[i]);
                        }
                    }
                }
                Reset(targetRegister);
            }
        }
        controlled (controls, ...) {
            Controlled ControlledXOnTruthTable (controls, (func, controlRegister, targetRegister));
        }
        controlled adjoint (controls, ...) {
            Controlled ControlledXOnTruthTable (controls, (func, controlRegister, targetRegister));
        }
    }
}
