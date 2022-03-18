// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Tests {
    open Microsoft.Quantum.Measurement;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Diagnostics as Diag;

    operation CheckAllowNCallsTestShouldFail() : Unit {
        within {
            Diag.AllowAtMostNCallsCA(3, X, "Too many calls to X.");
        } apply {
            use q = Qubit();
            // Should run four times, one more
            // than the three allowed.
            for idx in 0..3 {
                X(q);
            }
        }
    }

    @Diag.Test("QuantumSimulator")
    operation CheckAllowNCalls() : Unit {
        within {
            Diag.AllowAtMostNCallsCA(4, X, "Too many calls to X.");
        } apply {
            use q = Qubit();
            // Should run four times, exactly as
            // many times as allowed.
            for idx in 0..3 {
                X(q);
            }
        }
    }

    operation CheckAllowNCallsMeasurementsTestShouldFail() : Unit {
        within {
            Diag.AllowAtMostNCallsCA(3, Measure, "Too many calls to Measure.");
        } apply {
            use q = Qubit();
            // Should use four measurements, one more
            // than the three allowed.
            // M should be expressed via Measure, so all calls should add up.
            let resM = M(q);
            let resMResetX = MResetX(q);
            let resMResetY = MResetY(q);
            let resMeasure = Measure([PauliZ], [q]);
        }
    }

    @Diag.Test("QuantumSimulator")
    operation CheckAllowNCallsMeasurements() : Unit {
        within {
            Diag.AllowAtMostNCallsCA(4, Measure, "Too many calls to Measure.");
        } apply {
            use q = Qubit();
            // Should use four measurements, exactly as
            // many times as allowed.
            // M should be expressed via Measure, so all calls should add up.
            let resM = M(q);
            let resMResetX = MResetX(q);
            let resMResetY = MResetY(q);
            let resMeasure = Measure([PauliZ], [q]);
        }
    }

    @Diag.Test("QuantumSimulator")
    operation CheckAllowNCallsDistinguishControlled() : Unit {
        within {
            Diag.AllowAtMostNCallsCA(1, Controlled X, "Too many calls to Controlled X.");
        } apply {
            use (q1, q2) = (Qubit(), Qubit());
            // Should use one Controlled X, exactly as many times as allowed.
            // X will not be accounted for.
            X(q2);
            Controlled X([q1], q2);
            X(q2);
        }
    }

    @Diag.Test("QuantumSimulator")
    operation CheckAllowNCallsDistinguishAdjoint() : Unit {
        within {
            Diag.AllowAtMostNCallsCA(4, S, "Too many calls to S.");
            Diag.AllowAtMostNCallsCA(2, Adjoint S, "Too many calls to Adjoint S.");
            Diag.AllowAtMostNCallsCA(2, Controlled S, "Too many calls to Controlled S.");
            Diag.AllowAtMostNCallsCA(1, Adjoint Controlled S, "Too many calls to Adjoint Controlled S.");
        } apply {
            use (q1, q2) = (Qubit(), Qubit());
            // Should use two Adjoint S (one via Adjoint Controlled S), but exactly
            // one Adjoint Controlled S.
            S(q1);
            Controlled S([q1], q2);
            Adjoint Controlled S([q2], q1);
            Adjoint S(q1);

            Reset(q1);
            Reset(q2);
        }
    }

    @Diag.Test("ToffoliSimulator")
    operation CheckAllowNQubitsWithNestedCalls() : Unit {
        // Here, the total number of allocated qubits exceeds our policy,
        // but the number of those qubits allocated inside the policy
        // condition is still OK such that this should pass.
        use outer = Qubit[4];
        within {
            Diag.AllowAtMostNQubits(5, "Too many additional qubit allocations.");
        } apply {
            use qs = Qubit[2] {
                use qs2 = Qubit[1] { }
                use qs2 = Qubit[2] { }
            }
        }
    }

    operation CheckAllowNQubitsTestShouldFail() : Unit {
        within {
            Diag.AllowAtMostNQubits(3, "Too many additional qubit allocations.");
        } apply {
            use qs = Qubit[2] {
                use qs2 = Qubit[1] { }
                use qs2 = Qubit[2] { }
            }
        }
    }

    @Diag.Test("QuantumSimulator")
    operation CheckAllowNQubits() : Unit {
        within {
            Diag.AllowAtMostNQubits(3, "Too many additional qubit allocations.");
        } apply {
            use qs = Qubit[3];
        }
    }

}
