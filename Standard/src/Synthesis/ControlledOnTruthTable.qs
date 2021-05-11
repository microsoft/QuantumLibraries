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
    /// ```qsharp
    /// FastHadamardTransformed([1, 1, 1, -1]); // [2, 2, 2, -2]
    /// ```
    ///
    /// # Reference
    /// Frank Yates: The design and analysis of factorial experiments, in:
    /// Technical Communication No. 35, Imperial Bureau of Soil Science,
    /// London (1937)
    internal function FastHadamardTransformed(func : Int[]) : Int[] {
        let bits = BitSizeI(Length(func) - 1);
        mutable res = func;
        for m in 0..bits - 1 {
            mutable s = 1 <<< m;
            for i in 0..(2 * s)..Length(func) - 1 {
                mutable k = i + s;
                for j in i..i + s - 1 {
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
    /// ```qsharp
    /// Extended([2, 2, 2, -2]); // [2, 2, 2, -2, -2, -2, -2, 2]
    /// ```
    internal function Extended(spectrum : Int[]) : Int[] {
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
    /// ```qsharp
    /// Encoded([false, false, false, true]); // [1, 1, 1, -1]
    /// ```
    internal function Encoded(table : Bool[]) : Int[] {
        return Mapped(RMEncoding, table);
    }

    /// # Summary
    /// Adjusts truth table from array of Booleans according to number of variables
    ///
    /// A new array is returned of length `2^numVars`, possibly
    /// requiring to extend `table`'s size with `false` entries
    /// or truncating it to `2^numVars` elements.
    ///
    /// # Input
    /// ## table
    /// Truth table as array of truth values
    /// ## numVars
    /// Number of variables
    ///
    /// # Output
    /// Size adjusted truth table
    internal function SizeAdjustedTruthTable(table : Bool[], numVars : Int) : Bool[] {
        let numEntries = 2^numVars;
        if (numEntries < Length(table)) {
            return table[...numEntries - 1];
        } elif (numEntries > Length(table)) {
            return Padded(-numEntries, false, table);
        } else {
            return table;
        }
    }

    /// # Summary
    /// Applies the @"microsoft.quantum.intrinsic.x" operation on `target`, if the Boolean function `func` evaluates
    /// to true for the classical assignment in `controlRegister`.
    ///
    /// # Description
    /// The operation implements the unitary operation
    /// \begin{align}
    ///    U\ket{x}\ket{y} = \ket{x}\ket{y \oplus f(x)}
    /// \end{align}
    /// where $x$ and $y$ represent `controlRegister` and `target`, respectively.
    ///
    /// The Boolean function $f$ is represented as a truth table in terms of a big integer.
    /// For example, the majority function on three inputs is represented by the bitstring
    /// `11101000`, where the most significant bit `1` corresponds to the input assignment `(1, 1, 1)`,
    /// and the least significant bit `0` corresponds to the input assignment `(0, 0, 0)`.
    /// It can be represented by the big integer `0xE8L` in hexadecimal notation or as `232L`
    /// in decimal notation.  The `L` suffix indicates that the constant is of type `BigInt`.
    /// More details on this representation can also be found in the [truth tables kata](https://github.com/microsoft/QuantumKatas/tree/main/TruthTables).
    ///
    /// The implementation makes use of @"microsoft.quantum.intrinsic.cnot"
    /// and @"microsoft.quantum.intrinsic.r1" gates.
    ///
    /// # Input
    /// ## func
    /// Boolean truth table represented as big integer
    /// ## controlRegister
    /// Register of control qubits
    /// ## target
    /// Target qubit
    ///
    /// # See Also
    /// - Microsoft.Quantum.Synthesis.ApplyXControlledOnTruthTableWithCleanTarget
    ///
    /// # References
    /// - [*N. Schuch*, *J. Siewert*, PRL 91, no. 027902, 2003, arXiv:quant-ph/0303063](https://arxiv.org/abs/quant-ph/0303063)
    /// - [*Mathias Soeken*, *Martin Roetteler*, arXiv:2005.12310](https://arxiv.org/abs/2005.12310)
    operation ApplyXControlledOnTruthTable (func : BigInt, controlRegister : Qubit[], target : Qubit) : Unit {
        body (...) {
            let vars = Length(controlRegister);

            let maxValue = 1L <<< 2^vars;
            Fact(func >= 0L and func < maxValue, $"Argument func must be value from 0 to {maxValue}");

            let tt = BigIntAsBoolArray(func);
            let table = Encoded(SizeAdjustedTruthTable(BigIntAsBoolArray(func), vars));
            let spectrum = Extended(FastHadamardTransformed(table));

            let qubits = controlRegister + [target];

            HY(target);

            for i in 0..vars {
                let start = 1 <<< i;
                for (offset, ctrl) in GrayCode(i) {
                    R1Frac(spectrum[start + offset], vars + 1, qubits[i]);
                    if (i != 0) {
                        CNOT(qubits[ctrl], qubits[i]);
                    }
                }
            }

            H(target);
        }
        adjoint self;
        controlled (controls, ...) {
            use q = Qubit();
            within {
                ApplyXControlledOnTruthTableWithCleanTarget(func, controlRegister, q);
            } apply {
                Controlled X(controls + [q], target);
            }
        }
        controlled adjoint self;
    }

    /// # Summary
    /// Applies the @"microsoft.quantum.intrinsic.x" operation on `target`, if the Boolean function `func` evaluates
    /// to true for the classical assignment in `controlRegister`.
    ///
    /// # Description
    /// This operation implements a special case of @"microsoft.quantum.synthesis.applyxcontrolledontruthtable",
    /// in which the target qubit is known to be in the $\ket{0}$ state.
    ///
    /// The implementation makes use of @"microsoft.quantum.intrinsic.cnot"
    /// and @"microsoft.quantum.intrinsic.r1" gates.  The implementation of the
    /// adjoint operation is optimized and uses measurement-based uncomputation.
    ///
    /// # Input
    /// ## func
    /// Boolean truth table represented as big integer
    /// ## controlRegister
    /// Register of control qubits
    /// ## target
    /// Target qubit (must be in $\ket{0}$ state)
    ///
    /// # See Also
    /// - Microsoft.Quantum.Synthesis.ApplyXControlledOnTruthTable
    operation ApplyXControlledOnTruthTableWithCleanTarget (func : BigInt, controlRegister : Qubit[], target : Qubit) : Unit {
        body (...) {
            let vars = Length(controlRegister);

            let maxValue = PowL(2L, 2^vars);
            Fact(func >= 0L and func < maxValue, $"Argument func must be value from 0 to {maxValue}");
            AssertQubit(Zero, target);

            let tt = BigIntAsBoolArray(func);
            let table = Encoded(SizeAdjustedTruthTable(BigIntAsBoolArray(func), vars));
            let spectrum = FastHadamardTransformed(table);

            HY(target);

            for (offset, ctrl) in GrayCode(vars) {
                R1Frac(-spectrum[offset], vars + 1, target);
                CNOT(controlRegister[ctrl], target);
            }

            H(target);
        }
        adjoint (...) {
            let vars = Length(controlRegister);

            let maxValue = PowL(2L, 2^vars);
            Fact(func >= 0L and func < maxValue, $"Argument func must be value from 0 to {maxValue}");

            let tt = BigIntAsBoolArray(func);
            let table = Encoded(SizeAdjustedTruthTable(BigIntAsBoolArray(func), vars));
            let spectrum = FastHadamardTransformed(table);

            H(target);
            AssertMeasurementProbability([PauliZ], [target], One, 0.5, "Probability of the measurement must be 0.5", 1e-10);

            if IsResultOne(M(target)) {
                for i in 0..vars - 1 {
                    let start = 1 <<< i;
                    for (offset, ctrl) in GrayCode(i) {
                        R1Frac(spectrum[start + offset], vars, controlRegister[i]);
                        if (i != 0) {
                            CNOT(controlRegister[ctrl], controlRegister[i]);
                        }
                    }
                }
                Reset(target);
            }
        }
        controlled (controls, ...) {
            Controlled ApplyXControlledOnTruthTable (controls, (func, controlRegister, target));
        }
        controlled adjoint (controls, ...) {
            Controlled ApplyXControlledOnTruthTable (controls, (func, controlRegister, target));
        }
    }
}
