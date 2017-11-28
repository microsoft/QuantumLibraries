// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {

    // TODO: what does C♯/LINQ call the (0, arr[0]), (1, arr[1]), ... IEnumerator?
    //       we might should use same naming convention here.
    operation ApplyToRangeAC(singleQubitOperation : ((Int, Qubit) => ():Adjoint,Controlled), register : Qubit[])  : ()
    {
        body {
            for (idxQubit in 0..Length(register) - 1) {
                singleQubitOperation(idxQubit, register[idxQubit]);
            }
        }

        adjoint auto
        controlled auto
        controlled adjoint auto
    }

}
