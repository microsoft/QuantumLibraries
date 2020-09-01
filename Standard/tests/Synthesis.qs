// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Tests {
    open Microsoft.Quantum.Arithmetic;
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Convert;
    open Microsoft.Quantum.Diagnostics;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Math;
    open Microsoft.Quantum.Measurement;
    open Microsoft.Quantum.Synthesis;
    open Microsoft.Quantum.Random;

    internal operation CheckApplyPermutation (synthesisOperation : ((Int[], LittleEndian) => Unit is Adj + Ctl)) : Unit {
        let permutations = [
            [0, 2, 1, 3],
            [0, 1, 3, 2],
            [0, 1, 2, 3],
            [3, 2, 1, 0],
            [0, 1, 2, 3, 4, 5, 7, 6],
            [0, 2, 4, 6, 1, 3, 5, 7],
            [0, 2, 3, 5, 7, 11, 13, 1, 4, 6, 8, 9, 10, 12, 14, 15]
        ];

        for (perm in permutations) {
            let numQubits = BitSizeI(Length(perm) - 1);

            using (qs = Qubit[numQubits]) {
                let register = LittleEndian(qs);
                for (i in 0..Length(perm) - 1) {
                    ApplyXorInPlace(i, register);
                    synthesisOperation(perm, register);
                    EqualityFactI(MeasureInteger(register), perm[i], $"ApplyPermutation failed for permutation {perm} at index {i}");
                }
            }
        }
    }

    internal operation ApplyPermutationUsingDecompositionWithReverseVariableOrder (perm : Int[], qubits : LittleEndian) : Unit is Adj + Ctl {
        let variableOrder = Reversed(SequenceI(0, Length(qubits!) - 1));
        ApplyPermutationUsingDecompositionWithVariableOrder(perm, variableOrder, qubits);
    }

    @Test("ToffoliSimulator")
    operation CheckTransformationBasedSynthesis () : Unit {
        CheckApplyPermutation(ApplyPermutationUsingTransformation);
    }

    @Test("QuantumSimulator")
    operation CheckDecompositionBasedSynthesis () : Unit {
        CheckApplyPermutation(ApplyPermutationUsingDecomposition);
        CheckApplyPermutation(ApplyPermutationUsingDecompositionWithReverseVariableOrder);
    }

    @Test("ToffoliSimulator")
    operation CheckApplyTransposition () : Unit {
        for (numQubits in 2..6) {
            for (_ in 1..10) {
                let a = DrawRandomInt(0, 2^numQubits - 1);
                let b = DrawRandomInt(0, 2^numQubits - 1);

                using (qs = Qubit[numQubits]) {
                    let register = LittleEndian(qs);

                    for (i in 0..2^numQubits - 1) {
                        ApplyXorInPlace(i, register);
                        ApplyTransposition(a, b, register);
                        EqualityFactI(
                            MeasureInteger(register),
                            i == a ? b | (i == b ? a | i),
                            $"ApplyTransposition failed for {numQubits} qubits when a = {a} and b = {b}"
                        );
                    }
                }
            }
        }
    }

    internal operation RandomBool() : Bool {
        return DrawRandomInt(0, 1) == 1;
    }

    // from ControlledOnTruthTable.qs
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

    @Test("QuantumSimulator")
    operation CheckApplyXControlledOnTruthTable () : Unit {
        for (numQubits in 2..5) {
            for (round in 1..5) {
                // TODO: replace with BigUniform.
                let func = IntAsBigInt(DrawRandomInt(0, 2^(2^numQubits) - 1));
                let truthValues = SizeAdjustedTruthTable(BigIntAsBoolArray(func), numQubits);

                using ((controls, target) = (Qubit[numQubits], Qubit())) {
                    for (i in 0..(2^numQubits - 1)) {
                        let targetInit = RandomBool();
                        ApplyIf(X, targetInit, target);
                        within {
                            ApplyXorInPlace(i, LittleEndian(controls));
                        } apply {
                            ApplyXControlledOnTruthTable(func, controls, target);
                        }
                        EqualityFactB(
                            IsResultOne(MResetZ(target)) != targetInit,
                            truthValues[i],
                            $"Measured value does not correspond to truth table bit in truth table {func} and bit {i}");
                    }
                }
            }
        }
    }

    @Test("QuantumSimulator")
    operation CheckControlledApplyXControlledOnTruthTable () : Unit {
        for (numQubits in 2..5) {
            for (round in 1..5) {
                // TODO: replace with BigUniform.
                let func = IntAsBigInt(DrawRandomInt(0, 2^(2^numQubits) - 1));
                let truthValues = SizeAdjustedTruthTable(BigIntAsBoolArray(func), numQubits);

                using ((controls, control, target) = (Qubit[numQubits], Qubit(), Qubit())) {
                    for (i in 0..(2^numQubits - 1)) {
                        let controlInit = RandomBool();
                        let targetInit = RandomBool();
                        ApplyIf(X, targetInit, target);
                        within {
                            ApplyIfA(X, controlInit, control);
                            ApplyXorInPlace(i, LittleEndian(controls));
                        } apply {
                            Controlled ApplyXControlledOnTruthTable([control], (func, controls, target));
                        }

                        let result = IsResultOne(MResetZ(target));
                        if (controlInit) {
                            EqualityFactB(
                                result != targetInit,
                                truthValues[i],
                                $"Measured value does not correspond to truth table bit in truth table {func} and bit {i}");
                        } else {
                            EqualityFactB(result, targetInit, $"Target should not have been changed from its initial value {targetInit}");
                        }
                    }
                }
            }
        }
    }

    @Test("QuantumSimulator")
    operation CheckApplyXControlledOnTruthTableWithCleanTarget () : Unit {
        for (numQubits in 2..5) {
            for (round in 1..5) {
                // TODO: replace with BigUniform.
                let func = IntAsBigInt(DrawRandomInt(0, 2^(2^numQubits) - 1));
                let truthValues = SizeAdjustedTruthTable(BigIntAsBoolArray(func), numQubits);

                using ((controls, target, copy) = (Qubit[numQubits], Qubit(), Qubit())) {
                    for (i in 0..(2^numQubits - 1)) {
                        within {
                            ApplyXorInPlace(i, LittleEndian(controls));
                            ApplyXControlledOnTruthTableWithCleanTarget(func, controls, target);
                        } apply {
                            CNOT(target, copy);
                        }
                        EqualityFactB(
                            IsResultOne(MResetZ(copy)),
                            truthValues[i],
                            $"Measured value does not correspond to truth table bit in truth table {func} and bit {i}");
                    }
                }
            }
        }
    }
}
