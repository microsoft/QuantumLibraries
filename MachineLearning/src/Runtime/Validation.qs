namespace Microsoft.Quantum.MachineLearning {
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Canon;

    function NMismatches(proposed: Int[], actual: Int[]): Int {
        mutable count = 0;
        for ((proposedLabel, actualLabel) in Zip(proposed, actual)) {
            if (proposedLabel != actualLabel) {
                set count += 1;
            }
        }
        return count;
    }

    /// # Summary
    /// tallies hits and misses off a list of probability estimates
    ///
    /// # Input
    /// ## pls
    /// a list of estimated probabilities with the corresponding class labels
    ///
    /// ## bias
    /// bias on record
    ///
    /// # Output
    /// (no.hits, no.misses) pair
    ///
    function TallyHitsMisses(pls : (Double, Int)[], bias : Double) : (Int, Int) {
        mutable hits = 0;
        mutable misses = 0;
        for ((classificationProbability, label) in pls) {
            if (label == InferredLabel(bias, classificationProbability)) {
                set hits += 1;
            } else {
                set misses += 1;
            }
        }
        return (hits, misses);
    }

    /// # Summary
    /// Using a flat description of a trained classification model, count
    /// the number of mispredictions occuring over the validation set
    ///
    /// # Input
    /// ## nQubits
    /// the number of qubits used for data encoding
    ///
    /// ## trainingSet
    /// the set of training samples
    ///
    /// ## trainingLabels
    /// the set of training labels
    ///
    /// ## validatioSchedule
    /// defines a subset of training data used for validation and computation of the *bias*
    ///
    /// ## gates
    /// Flat representation of classifier structure. Each element is
    /// [parameterIndex, pauliCode, targetQubit, sequence of control qubits]
    ///
    /// ## parameters
    /// an array of candidate parameters
    ///
    /// ## bias
    /// candidate predition bias
    ///
    /// ## nMeasurenets
    /// number of the measurement cycles to be used for estimation of each probability
    ///
    /// # Output
    /// the number of misclassifications
    ///
    operation CountValidationMisses(tolerance: Double, nQubits: Int, trainingSet: Double[][], trainingLabels: Int[], validationSchedule: Int[][], gates: Int[][], parameters: Double[],bias:Double, nMeasurements: Int) : Int
    {
        let schValidate = unFlattenSchedule(validationSchedule);
        let results = ValidateModel(
            tolerance, nQubits, Mapped(LabeledSample, Zip(trainingSet, trainingLabels)),
            schValidate, unFlattenGateSequence(gates),
            parameters, bias, nMeasurements
        );
        return results::NMisclassifications;
    }

    operation ValidateModel(tolerance: Double, nQubits: Int, samples : LabeledSample[], validationSchedule: SamplingSchedule, gates: GateSequence, parameters: Double[], bias:Double, nMeasurements: Int) : ValidationResults
    {
        let features = Mapped(_Features, samples);
        let labels = Sampled(validationSchedule, Mapped(_Label, samples));
        let probsValidation = EstimateClassificationProbabilitiesClassicalData(tolerance, features, validationSchedule, nQubits,  gates, parameters, nMeasurements);
        let localPL = InferredLabels(bias, probsValidation);
        let nMismatches = NMismatches(localPL, labels);
        return ValidationResults(
            nMismatches
        );
    }

}
