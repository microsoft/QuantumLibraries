// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Tests {
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Diagnostics as Diag;

    operation DumpS() : Unit {
        Diag.DumpOperation(1, ApplyToEachCA(S, _));
    }

    operation ApplyCnotToRegister(register : Qubit[]) : Unit is Adj + Ctl {
        CNOT(register[0], register[1]);
    }

    operation DumpCnot() : Unit {
        Diag.DumpOperation(2, ApplyCnotToRegister);
    }

}
