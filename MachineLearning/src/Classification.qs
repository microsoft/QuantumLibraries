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
        model : SequentialModel,
        target : Qubit[]
    )
    : Unit is Adj {
        encoder(LittleEndian(target));
        ApplySequentialClassifier(model, target);
    }

    /// # Summary
    /// Given a sample and a sequential classifier, estimates the
    /// classification probability for that sample by repeatedly measuring
    /// the output of the classifier on the given sample.
    ///
    /// # Input
    /// ## tolerance
    /// The tolerance to allow in encoding the sample into a state preparation
    /// operation.
    /// ## parameters
    /// A parameterization of the given sequential classifier.
    /// ## structure
    /// The structure of the given sequential classifier.
    /// ## sample
    /// The feature vector for the sample to be classified.
    /// ## nMeasurements
    /// The number of measusrements to use in estimating the classification
    /// probability.
    /// # Output
    /// An estimate of the classification probability for the given sample.
    operation EstimateClassificationProbability(
        tolerance : Double,
        model : SequentialModel,
        sample : Double[],
        nMeasurements: Int
    )
    : Double {
        let nQubits = FeatureRegisterSize(sample);
        let circEnc = ApproximateInputEncoder(tolerance / IntAsDouble(Length(model::Structure)), sample);
        let encodedSample = StateGenerator(nQubits, circEnc);
        return 1.0 - EstimateFrequencyA(
            _PrepareClassification(encodedSample::Prepare, model, _),
            _TailMeasurement(encodedSample::NQubits),
            encodedSample::NQubits,
            nMeasurements
        );
    }

    /// # Summary
    /// Given a set of samples and a sequential classifier, estimates the
    /// classification probability for those samples by repeatedly measuring
    /// the output of the classifier on each sample.
    ///
    /// # Input
    /// ## tolerance
    /// The tolerance to allow in encoding the sample into a state preparation
    /// operation.
    /// ## parameters
    /// A parameterization of the given sequential classifier.
    /// ## structure
    /// The structure of the given sequential classifier.
    /// ## samples
    /// An array of feature vectors for each sample to be classified.
    /// ## nMeasurements
    /// The number of measusrements to use in estimating the classification
    /// probability.
    /// # Output
    /// An array of estimates of the classification probability for each given
    /// sample.
    operation EstimateClassificationProbabilities(
        tolerance : Double,
        model : SequentialModel,
        samples : Double[][],
        nMeasurements : Int
    )
    : Double[] {
        let effectiveTolerance = tolerance / IntAsDouble(Length(model::Structure));
        return ForEach(
            EstimateClassificationProbability(
                effectiveTolerance, model, _, nMeasurements
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
