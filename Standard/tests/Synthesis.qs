// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.
namespace Microsoft.Quantum.Tests {
    open Microsoft.Quantum.Arithmetic;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Convert;
    open Microsoft.Quantum.Diagnostics;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Math;
    open Microsoft.Quantum.Measurement;
    open Microsoft.Quantum.Synthesis;

    @Test("ToffoliSimulator")
    operation CheckTransformationBasedSynthesis () : Unit {
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
                    ApplyPermutationUsingTransformation(perm, register);
                    EqualityFactI(MeasureInteger(register), perm[i], $"ApplyPermutation failed for permutation {perm} at index {i}");
                }
            }
        }
    }

    @Test("QuantumSimulator")
    operation CheckControlledXOnTruthTable () : Unit {
        for (numQubits in 2..5) {
            for (round in 1..5) {
                let func = IntAsBigInt(RandomInt(2^(2^numQubits)));
                let truthValues = BigIntAsBoolArray(func);

                using ((controls, target) = (Qubit[numQubits], Qubit())) {
                    for (i in 0..(2^numQubits - 1)) {
                        within {
                            ApplyXorInPlace(i, LittleEndian(controls));
                        } apply {
                            ControlledXOnTruthTable(func, controls, target);
                        }
                        EqualityFactB(
                            IsResultOne(MResetZ(target)),
                            truthValues[i],
                            $"Measured value does not correspond to truth table bit in truth table {func} and bit {i}");
                    }
                }
            }
        }
    }
}
