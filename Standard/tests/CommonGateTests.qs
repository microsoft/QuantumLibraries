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
        using (ancilla = Qubit[2]) {
            // 00
            ApplyFermionicSWAP(ancilla[0], ancilla[1]);
            AssertAllZero(ancilla);

            // 01
            X(ancilla[1]);
            ApplyFermionicSWAP(ancilla[0], ancilla[1]);
            X(ancilla[0]);
            AssertAllZero(ancilla);

            // 10
            X(ancilla[0]);
            ApplyFermionicSWAP(ancilla[0], ancilla[1]);
            X(ancilla[1]);
            AssertAllZero(ancilla);

            // 11
            ApplyToEachCA(X, ancilla);
            ApplyFermionicSWAP(ancilla[0], ancilla[1]);
            ApplyToEachCA(X, ancilla);
            AssertAllZero(ancilla);
        }
    }

    operation VerifyFermionicSWAPPhaseHelper(qubit1 : Qubit, qubit2: Qubit, phase : Result) : Unit {
        ApplyFermionicSWAP(qubit1, qubit2);
        Assert([PauliZ, PauliZ], [qubit1, qubit2], phase,
            "The Fermionic SWAP applies an incorrect phase");
    }
    
    // Verify Fermionic SWAP gives the correct phase change
    operation ApplyFermionicSWAPPhaseTest() : Unit {
        using (ancilla = Qubit[2]) {
            // 00
            VerifyFermionicSWAPPhaseHelper(ancilla[0], ancilla[1], Zero);
            ResetAll(ancilla);

            // 01
            X(ancilla[1]);
            VerifyFermionicSWAPPhaseHelper(ancilla[0], ancilla[1], One);
            ResetAll(ancilla);

            // 10
            X(ancilla[0]);
            VerifyFermionicSWAPPhaseHelper(ancilla[0], ancilla[1], One);
            ResetAll(ancilla);

            // 11
            ApplyToEachCA(X, ancilla);
            VerifyFermionicSWAPPhaseHelper(ancilla[0], ancilla[1], Zero);
            ResetAll(ancilla);
        }
    }

}


