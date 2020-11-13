// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.AmplitudeAmplification {

    /// # Summary
    /// Phases for a sequence of partial reflections in amplitude amplification.
    ///
    /// # Named Items
    /// ## AboutStart
    /// An array of phases for reflection about the
    /// start state.
    /// ## AboutTarget
    /// An array of phases for reflection
    /// about the target state.
    ///
    /// # Remarks
    /// Both arrays must be of equal length. Note that in many cases, the first phase about the start state and last phase about the target state introduces a global phase shift and may be set to $0$.
    newtype ReflectionPhases = (
        AboutStart: Double[],
        AboutTarget: Double[]
    );

    /// # Summary
    /// Phases for a sequence of single-qubit rotations in amplitude amplification.
    ///
    /// # Remarks
    /// The first parameter is an array of phases for reflections, expressed as a product of single-qubit rotations.
    /// [ G.H. Low, I. L. Chuang, https://arxiv.org/abs/1707.05391].
    newtype RotationPhases = Double[];

}


