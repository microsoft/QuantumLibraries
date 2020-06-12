// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.
namespace Microsoft.Quantum.Tests {
    open Microsoft.Quantum.Arithmetic;
    open Microsoft.Quantum.Canon;
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
            [0, 2, 3, 5, 7, 11, 13, 4, 6, 8, 9, 10, 12, 14, 15]
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
}
