// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Tests {
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Diagnostics as Diag;

    operation CheckAllowNCallsTestShouldFail() : Unit {
        within {
            Diag.AllowAtMostNCallsCA(3, X);
        } apply {
            using (q = Qubit()) {
                // Should run four times, one more
                // than the three allowed.
                for (idx in 0..3) {
                    X(q);
                }
            }
        }
    }

    @Diag.Test("QuantumSimulator")
    operation CheckAllowNCalls() : Unit {
        within {
            Diag.AllowAtMostNCallsCA(4, X);
        } apply {
            using (q = Qubit()) {
                // Should run four times, exactly as
                // many times as allowed.
                for (idx in 0..3) {
                    X(q);
                }
            }
        }
    }

    operation CheckAllowNQubitsTestShouldFail() : Unit {
        within {
            Diag.AllowAtMostNQubits(3);
        } apply {
            using (qs = Qubit[2]) {
                using (qs2 = Qubit[1]) { }
                using (qs2 = Qubit[2]) { }
            }
        }
    }

    @Diag.Test("QuantumSimulator")
    operation CheckAllowNQubits() : Unit {
        within {
            Diag.AllowAtMostNQubits(3);
        } apply {
            using (qs = Qubit[3]) { }
        }
    }

}
