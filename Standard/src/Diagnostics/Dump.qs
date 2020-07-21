// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Diagnostics {
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Characterization;
    open Microsoft.Quantum.Preparation;

    operation DumpOperation(nQubits : Int, op : (Qubit[] => Unit)) : Unit {
        using ((reference, target) = (Qubit[nQubits], Qubit[nQubits])) {
            PrepareChoiState(op, reference, target);
            DumpReferenceAndTarget(reference, target);
            ResetAll(reference + target);
        }
    }

}
