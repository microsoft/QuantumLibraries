// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.AmplitudeAmplification {
    open Microsoft.Quantum.Math;

    /// # Summary
    /// Converts phases specified as single-qubit rotations to phases
    /// specified as partial reflections.
    ///
    /// # Input
    /// ## rotPhases
    /// Array of single-qubit rotations to be converted to partial
    /// reflections.
    ///
    /// # Output
    /// An operation that implements phases specified as partial reflections.
    ///
    /// # References
    /// We use the convention in
    /// - [ *G.H. Low, I. L. Chuang* ](https://arxiv.org/abs/1707.05391)
    /// for relating single-qubit rotation phases to reflection operator phases.
    function RotationPhasesAsReflectionPhases(rotPhases : RotationPhases)
    : ReflectionPhases {
        let nPhasesRot = Length(rotPhases!);
        let nPhasesRef = (nPhasesRot + 1) / 2;

        if (nPhasesRot % 2 == 0) {
            fail $"Number of rotations must be odd.";
        }

        mutable phasesTarget = [0.0, size = nPhasesRef];
        mutable phasesStart = [0.0, size = nPhasesRef];

        if nPhasesRot == 1 {
            return ReflectionPhases(phasesStart, phasesTarget);
        }

        set phasesTarget w/= 0 <- ((rotPhases!)[0] - (rotPhases!)[1]) - PI();
        set phasesStart w/= 0 <- -(rotPhases!)[0] + 0.5 * PI();

        for idxPhases in 1 .. nPhasesRef - 2 {
            set phasesTarget w/= idxPhases <- ((rotPhases!)[2 * idxPhases] - (rotPhases!)[2 * idxPhases + 1]) - PI();
            set phasesStart w/= idxPhases <- ((rotPhases!)[2 * idxPhases - 1] - (rotPhases!)[2 * idxPhases]) + PI();
        }

        set phasesTarget w/= nPhasesRef - 1 <- (rotPhases!)[2 * nPhasesRef - 2] - 0.5 * PI();
        set phasesStart w/= nPhasesRef - 1 <- ((rotPhases!)[2 * nPhasesRef - 3] - (rotPhases!)[2 * nPhasesRef - 2]) + PI();
        return ReflectionPhases(phasesStart, phasesTarget);
    }

}
