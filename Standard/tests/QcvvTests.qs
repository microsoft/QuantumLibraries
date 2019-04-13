// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.
namespace Microsoft.Quantum.Tests {
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Extensions.Convert;
    open Microsoft.Quantum.Oracles;
    open Microsoft.Quantum.Characterization;
    open Microsoft.Quantum.Preparation;
    open Microsoft.Quantum.Diagnostics;
    open Microsoft.Quantum.Arrays;

    operation ChoiStateTest () : Unit {
        using (register = Qubit[2]) {
            PrepareChoiStateCA(NoOp<Qubit[]>, [register[0]], [register[1]]);

            // As usual, the same confusion about {+1, -1} and {0, 1}
            // labeling bites us here.
            Assert([PauliX, PauliX], register, Zero, $"XX");
            Assert([PauliZ, PauliZ], register, Zero, $"ZZ");
            ResetAll(register);
        }
    }

    operation EstimateFrequencyTest () : Unit {
        let freq = EstimateFrequency(ApplyToEach(H, _), MeasureAllZ, 1, 1000);
        EqualityWithinToleranceFact(freq, 0.5, 0.1);
    }

    operation _RobustPhaseEstimationTestOp (phase : Double, power : Int, qubits : Qubit[]) : Unit is Adj + Ctl {
        Exp([PauliZ], phase * ToDouble(power), qubits);
    }

    operation RobustPhaseEstimationDemoImpl (phaseSet : Double, bitsPrecision : Int) : Double {
        let op = DiscreteOracle(_RobustPhaseEstimationTestOp(phaseSet, _, _));

        using (q = Qubit()) {
            let phaseEst = RobustPhaseEstimation(bitsPrecision, op, [q]);
            Reset(q);
            return phaseEst;
        }
    }

    // Probabilistic test. Might fail occasionally
    operation RobustPhaseEstimationTest () : Unit {
        
        let bitsPrecision = 10;
        
        for (idxTest in 0 .. 9) {
            let phaseSet = ((2.0 * Microsoft.Quantum.Extensions.Math.PI()) * ToDouble(idxTest - 5)) / 12.0;
            let phaseEst = RobustPhaseEstimationDemoImpl(phaseSet, bitsPrecision);
            EqualityWithinToleranceFact(phaseEst, phaseSet, 0.01);
        }
    }
    
    
    operation PrepareQubitTest () : Unit {
        using (qubit = Qubit()) {
            let bases = [PauliI, PauliX, PauliY, PauliZ];

            for (basis in bases) {
                PrepareQubit(basis, qubit);
                Assert([basis], [qubit], Zero, $"Did not prepare in {basis} correctly.");
                Reset(qubit);
            }
        }
    }

    operation SingleQubitProcessTomographyMeasurementTest () : Unit {
        EqualityFactR(SingleQubitProcessTomographyMeasurement(PauliI, PauliI, H), Zero, $"Failed at ⟪I | H | I⟫.");
        EqualityFactR(SingleQubitProcessTomographyMeasurement(PauliX, PauliI, H), Zero, $"Failed at ⟪I | H | X⟫.");
        EqualityFactR(SingleQubitProcessTomographyMeasurement(PauliY, PauliI, H), Zero, $"Failed at ⟪I | H | Y⟫.");
        EqualityFactR(SingleQubitProcessTomographyMeasurement(PauliZ, PauliI, H), Zero, $"Failed at ⟪I | H | Z⟫.");
        EqualityFactR(SingleQubitProcessTomographyMeasurement(PauliX, PauliZ, H), Zero, $"Failed at ⟪Z | H | X⟫.");
        EqualityFactR(SingleQubitProcessTomographyMeasurement(PauliY, PauliY, H), One, $"Failed at -⟪Y | H | Y⟫.");
        EqualityFactR(SingleQubitProcessTomographyMeasurement(PauliX, PauliZ, H), Zero, $"Failed at ⟪Z | H | X⟫.");
    }

}


