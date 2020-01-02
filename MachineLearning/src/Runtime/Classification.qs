namespace Microsoft.Quantum.MachineLearning {
    open Microsoft.Quantum.Math;
    open Microsoft.Quantum.Characterization;
    open Microsoft.Quantum.Arithmetic;
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Convert;

    
    operation _PrepareClassification(
        encoder : (LittleEndian => Unit is Adj + Ctl),
        parameters : Double[],
        gates : GateSequence,
        target : Qubit[]
    )
    : Unit is Adj {
        encoder(LittleEndian(target));
        _ApplyGates(parameters, gates, target);
    }

    operation EstimateClassificationProbability(
        tolerance: Double,
        parameters : Double[],
        gates: GateSequence,
        sample: Double[],
        nMeasurements: Int
    )
    : Double {
        let nQubits = FeatureRegisterSize(sample);
        let circEnc = NoisyInputEncoder(tolerance / IntAsDouble(Length(gates!)), sample);
        let encodedSample = StateGenerator(nQubits, circEnc);
        return 1.0 - EstimateFrequencyA(
            _PrepareClassification(encodedSample::Apply, parameters, gates, _),
            _TailMeasurement(encodedSample::NQubits),
            encodedSample::NQubits,
            nMeasurements
        );
    }

    operation EstimateClassificationProbabilities(
        tolerance : Double,
        parameters : Double[],
        structure : GateSequence,
        samples : Double[][],
        nMeasurements : Int
    )
    : Double[] {
        let effectiveTolerance = tolerance / IntAsDouble(Length(structure!));
        return ForEach(
            EstimateClassificationProbability(
                effectiveTolerance, parameters, structure, _, nMeasurements
            ),
            samples
        );
    }

    /// # Summary
    /// Given a of classification probability and a bias, returns the
    /// label inferred from that probability.
    ///
    /// # Input
    /// ## bias
    /// The bias between two classes, typically the result of training a
    /// classifier.
    /// ## probability
    /// A classification probabilities for a particular sample, typicaly
    /// resulting from estimating its classification frequency.
    ///
    /// # Output
    /// The label inferred from the given classification probability.
    function InferredLabel(bias : Double, probability : Double) : Int {
        return probability + bias > 0.5 ? 1 | 0;
    }

    /// # Summary
    /// Given an array of classification probabilities and a bias, returns the
    /// label inferred from each probability.
    ///
    /// # Input
    /// ## bias
    /// The bias between two classes, typically the result of training a
    /// classifier.
    /// ## probabilities
    /// An array of classification probabilities for a set of samples, typicaly
    /// resulting from estimating classification frequencies.
    ///
    /// # Output
    /// The label inferred from each classification probability.
    function InferredLabels(bias : Double, probabilities : Double[]): Int[] {
        return Mapped(InferredLabel(bias, _), probabilities);
    }

}
