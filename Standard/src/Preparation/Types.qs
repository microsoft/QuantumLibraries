// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Preparation {
    open Microsoft.Quantum.Arithmetic;

    newtype MixedPreparationOperation = (
        Requirements: MixedPreparationRequirements,
        Norm: Double,
        Prepare: ((LittleEndian, Qubit[]) => Unit is Adj + Ctl)
    );

    newtype MixedPreparationRequirements = (
        NTotalQubits: Int,
        (
            NIndexQubits: Int,
            NGarbageQubits: Int
        )
    );

}
