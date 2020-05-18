// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.
namespace Microsoft.Quantum.Tests {
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Convert;
    open Microsoft.Quantum.Oracles;
    open Microsoft.Quantum.Characterization;
    open Microsoft.Quantum.Preparation;
    open Microsoft.Quantum.Diagnostics;
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Math;

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

    internal operation PrepareTrivialState(qubits : Qubit[]) : Unit is Adj {
        Message("stage prepared");
    }

    operation EstimateFrequencyTest () : Unit {
        let freq1 = EstimateFrequency(ApplyToEach(H, _), MeasureAllZ, 1, 1000);
        EqualityWithinToleranceFact(freq1, 0.5, 0.1);

        let freq2 = EstimateFrequencyA(ApplyToEachA(H, _), MeasureAllZ, 3, 10000);
        EqualityWithinToleranceFact(freq2, 0.5, 0.1);
    }

    operation PrepareBiasedCoin(successProbability : Double, qubit : Qubit) : Unit is Adj {
        let rotationAngle = 2.0 * ArcCos(Sqrt(successProbability));
        Ry(rotationAngle, qubit);
    }

    operation EstimateFrequencyBinomialCase(nSamples : Int, successProbability : Double, nStandardDeviations : Double) : Unit {
        let expectation = successProbability;
        let tolerance = nStandardDeviations * Sqrt(
            (successProbability * (1.0 - successProbability)) / IntAsDouble(nSamples)
        );
        let actualFreq = EstimateFrequencyA(
            ApplyToEachA(PrepareBiasedCoin(successProbability, _), _),
            Measure([PauliZ], _),
            1,
            nSamples
        );
        EqualityWithinToleranceFact(expectation, actualFreq, tolerance);
    }

    operation EstimateFrequencyBinomialTest() : Unit {
        // If this is larger, tests fail less often, but more false negatives
        // slip through.
        let nStdDevs = 3.0;
        for (testCase in [
            // 𝑛𝑝 <= 30
            (45, 0.5, nStdDevs),
            (100, 0.2, nStdDevs),
            (1000, 0.02, nStdDevs),
            // 𝑛𝑝 > 30
            (10000, 0.5, nStdDevs),
            (100000, 0.5, nStdDevs),
            (100000, 0.3, nStdDevs),
            (100000, 0.95, nStdDevs)
        ]) {
            EstimateFrequencyBinomialCase(testCase);
        }
    }

    // Calls EstimateFrequency with a TrivialStatePreparation to make sure
    // Emulation is actually kicking in.
    @Test("QuantumSimulator")
    @Test("ToffoliSimulator")
    operation EstimateFrequencyEmulationTest() : Unit {
        let freq = EstimateFrequencyA(PrepareTrivialState, Measure([PauliZ, PauliZ, PauliZ], _), 3, 2000);
        NearEqualityFactD(freq, 1.0);
    }

    internal operation RobustPhaseEstimationTestOp (phase : Double, power : Int, qubits : Qubit[]) : Unit is Adj + Ctl {
        Exp([PauliZ], phase * IntAsDouble(power), qubits);
    }

    operation RobustPhaseEstimationDemoImpl (phaseSet : Double, bitsPrecision : Int) : Double {
        let op = DiscreteOracle(RobustPhaseEstimationTestOp(phaseSet, _, _));

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
            let phaseSet = ((2.0 * PI()) * IntAsDouble(idxTest - 5)) / 12.0;
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
