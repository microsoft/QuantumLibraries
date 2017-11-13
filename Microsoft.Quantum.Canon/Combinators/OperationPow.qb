// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {

    operation OperationPowImpl(oracle : (Qubit[] => ()), power : Int, target : Qubit[])  : ()
    {
        body {
            for (idxApplication in 0..power - 1) {
                oracle(target);
            }
        }
    }

    operation OperationPowImplAC(oracle : (Qubit[] => ():Controlled,Adjoint), power : Int, target : Qubit[])  : ()
    {
        body {
            for (idxApplication in 0..power - 1) {
                oracle(target);
            }
        }

        adjoint auto
        controlled auto
        controlled adjoint auto
    }

    operation OperationPow(oracle : (Qubit[] => ()), power : Int)  : (Qubit[] => ())
    {
        body {
            return OperationPowImpl(oracle, power, _);
        }
    }

}
