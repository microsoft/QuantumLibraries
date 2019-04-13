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

    operation ApplyControlledOpToRegister(op : (Qubit => Unit : Adjoint, Controlled), target : Qubit[]) : Unit {
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
    
}


