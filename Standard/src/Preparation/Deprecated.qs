// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Preparation {

    @Deprecated("Microsoft.Quantum.Preparation.PrepareSingleQubitPositivePauliEigenstate")
    operation PrepareQubit (basis : Pauli, qubit : Qubit) : Unit {
        PrepareSingleQubitPositivePauliEigenstate(basis, qubit);
    }
}
