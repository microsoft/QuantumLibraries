// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.
namespace Microsoft.Quantum.Tests {
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Diagnostics;
    open Microsoft.Quantum.Arrays;

    operation ApplyCShorthandToRegister(cGate : ((Qubit, Qubit) => Unit), target : Qubit[]) : Unit {
        cGate(target[0], target[1]);
    }

    operation ApplyControlledOpToRegister(op : (Qubit => Unit is Adj + Ctl), target : Qubit[]) : Unit {
        body (...) {
            Controlled op(Most(target), Tail(target));
        }
        adjoint auto;
    }
    
    operation CXTest() : Unit {
        let actual = ApplyCShorthandToRegister(CX, _);
        let expected = ApplyControlledOpToRegister(X, _);
        AssertOperationsEqualReferenced(2, actual, expected);
    }

    operation CYTest() : Unit {
        let actual = ApplyCShorthandToRegister(CY, _);
        let expected = ApplyControlledOpToRegister(Y, _);
        AssertOperationsEqualReferenced(2, actual, expected);
    }

    operation CZTest() : Unit {
        let actual = ApplyCShorthandToRegister(CZ, _);
        let expected = ApplyControlledOpToRegister(Z, _);
        AssertOperationsEqualReferenced(2, actual, expected);
    }

    // Verify Fermionic SWAP gives the correct qubit values
    operation ApplyFermionicSWAPValueTest() : Unit {
        using ((left, right) = (Qubit(), Qubit())) {
            // 00
            ApplyFermionicSWAP(left, right);
            AssertAllZero([left, right]);

            // 01
            X(right);
            ApplyFermionicSWAP(left, right);
            X(left);
            AssertAllZero([left, right]);

            // 10
            X(left);
            ApplyFermionicSWAP(left, right);
            X(right);
            AssertAllZero([left, right]);

            // 11
            ApplyToEachCA(X, [left, right]);
            ApplyFermionicSWAP(left, right);
            ApplyToEachCA(X, [left, right]);
            AssertAllZero([left, right]);
        }
    }

    operation VerifyFermionicSWAPPhaseHelper(phase : Result, qubit1 : Qubit, qubit2: Qubit) : Unit {
        ApplyFermionicSWAP(qubit1, qubit2);
        AssertMeasurement([PauliZ, PauliZ], [qubit1, qubit2], phase,
            "The Fermionic SWAP applies an incorrect phase");
    }
    
    // Verify Fermionic SWAP gives the correct phase change
    operation ApplyFermionicSWAPPhaseTest() : Unit {
        using ((left, right) = (Qubit(), Qubit())) {
            // 00
            VerifyFermionicSWAPPhaseHelper(Zero, left, right);
            ResetAll([left, right]);

            // 01
            X(right);
            VerifyFermionicSWAPPhaseHelper(One, left, right);
            ResetAll([left, right]);

            // 10
            X(left);
            VerifyFermionicSWAPPhaseHelper(One, left, right);
            ResetAll([left, right]);

            // 11
            ApplyToEachCA(X, [left, right]);
            VerifyFermionicSWAPPhaseHelper(Zero, left, right);
            ResetAll([left, right]);
        }
    }

}


