// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {

    @Deprecated("Microsoft.Quantum.Canon.DecomposedIntoTimeStepsCA")
    function DecomposeIntoTimeStepsCA<'T>(
        (nSteps : Int, op : ((Int, Double, 'T) => Unit is Adj + Ctl)),
        trotterOrder : Int
    )
    : ((Double, 'T) => Unit is Adj + Ctl) {
        return DecomposedIntoTimeStepsCA((nSteps, op), trotterOrder);
    }

}


