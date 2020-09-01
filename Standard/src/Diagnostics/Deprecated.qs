// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {

    @Deprecated("Microsoft.Quantum.Diagnostics.AssertPhase")
    operation AssertPhase(expected : Double, qubit : Qubit, tolerance : Double) : Unit {
        Microsoft.Quantum.Diagnostics.AssertPhase(expected, qubit, tolerance);
    }

}
