namespace Microsoft.Quantum.Tests {
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Arithmetic;
    open Microsoft.Quantum.Convert;
    open Microsoft.Quantum.Oracles;
    open Microsoft.Quantum.Characterization;
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Math;
    open Microsoft.Quantum.Diagnostics;


    /// # Summary
    /// Assert that the QuantumPhaseEstimation operation for the T gate
    /// return 0000 in the controlRegister when targetState is 0 and
    /// return 0010 when the targetState is 1
    operation QuantumPhaseEstimationTest () : Unit {

        let oracle = DiscreteOracle(ApplyTOracle);

        using (qPhase = Qubit[5]) {
            let phase = BigEndian(qPhase[0 .. 3]);
            let state = qPhase[4];
            QuantumPhaseEstimation(oracle, [state], phase);
            let complexOne = Complex(1.0, 0.0);
            let complexZero = Complex(0.0, 0.0);

            for (idxPhase in 0 .. 4) {
                AssertQubitIsInStateWithinTolerance((complexOne, complexZero), qPhase[idxPhase], 1E-06);
            }

            X(state);
            QuantumPhaseEstimation(oracle, [state], phase);
            AssertQubitIsInStateWithinTolerance((complexOne, complexZero), qPhase[0], 1E-06);
            AssertQubitIsInStateWithinTolerance((complexOne, complexZero), qPhase[1], 1E-06);
            AssertQubitIsInStateWithinTolerance((complexZero, complexOne), qPhase[2], 1E-06);
            AssertQubitIsInStateWithinTolerance((complexOne, complexZero), qPhase[3], 1E-06);
            AssertQubitIsInStateWithinTolerance((complexZero, complexOne), qPhase[4], 1E-06);
            ResetAll(qPhase);
        }
    }


    /// # Summary
    /// Implementation of T-gate for Quantum Phase Estimation Oracle
    operation ApplyTOracle (power : Int, target : Qubit[]) : Unit is Adj + Ctl {
        for (idxPower in 0 .. power - 1) {
            T(Head(target));
        }
    }

}


