// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Tests {
    
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Diagnostics;
    open Microsoft.Quantum.Arrays;


    operation ApplySeriesOfOpsTest() : Unit {
        // create the sample ops + their targets here
        let op1 = ApplyToFirstQubit(X, _);
        let op2 = ApplyToFirstTwoQubits(CNOT, _);
        let op3 = Exp([PauliX, PauliZ, PauliY], 0.2002, _);
        let op4 = ApplyToEachA(H, _);
        let target1 = [0];
        let target2 = [0, 4];
        let target3 = [2, 3, 5];
        let target4 = [1, 2, 3, 4];

        let listOfOps = [op1, op2, op3, op4];
        let listOfTargets = [target1, target2, target3, target4];
        AssertOperationsEqualReferenced(6, ApplySeriesOfOps(listOfOps, listOfTargets, _), _SampleApplySeriesOfOps(_));
    }

    // Helper method for ApplySeriesOfOpsTest
    operation _SampleApplySeriesOfOps(register : Qubit[]) : Unit is Adj + Ctl {
        // replicate those ops implemented here
        X(register[0]);
        CNOT(register[0], register[4]);
        Exp([PauliX, PauliZ, PauliY], 0.2002, Subarray([2, 3, 5], register));
        ApplyToEachA(H, Subarray([1, 2, 3, 4], register));
    }
    
    operation ApplyRepeatedOpTest() : Unit {
        let op = ApplyToFirstThreeQubits(CCNOT, _);
        let targets = [[0, 1, 2], [2, 1, 0], [3, 4, 5], [2, 4, 0], [5, 3, 1]];
        AssertOperationsEqualReferenced(6, ApplyOpRepeatedlyOver(op, targets, _), _SampleApplyRepeatedOp(_));
    }

    // Helper method for ApplyRepeatedOpTest
    operation _SampleApplyRepeatedOp(register : Qubit[]) : Unit is Adj + Ctl {
        CCNOT(register[0], register[1], register[2]);
        CCNOT(register[2], register[1], register[0]);
        CCNOT(register[3], register[4], register[5]);
        CCNOT(register[2], register[4], register[0]);
        CCNOT(register[5], register[3], register[1]);
    }

    operation PermuteQubitsTest() : Unit {
        let sampleOrder = [5, 3, 2, 0, 1, 4];
        AssertOperationsEqualReferenced(6, PermuteQubits(sampleOrder, _) , _SamplePermuteQubits);
    }

    // Helper method for PermuteQubitsTest
    operation _SamplePermuteQubits(register : Qubit[]) : Unit is Adj + Ctl {
        // assumes the order to be swapped is [(0, 5),(0, 4),(0, 1),(0, 3)]
        // (Order is [5, 3, 2, 0, 1, 4])
        SWAP(register[0], register[5]);
        SWAP(register[0], register[4]);
        SWAP(register[0], register[1]);
        SWAP(register[0], register[3]);
    } 
    
}