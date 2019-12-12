namespace Microsoft.Quantum.MachineLearning {
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Canon;

    function MissLocations(inferredLabels : Int[], actualLabels : Int[])
    : Int[] {
        mutable ret = new Int[0];
        for ((idx, (inferred, actual)) in Enumerated(Zip(inferredLabels, actualLabels))) {
            if (inferred != actual) {
                set ret += [idx];
            }
        }
        return ret;
    }

    function NMismatches(proposed: Int[], actual: Int[]): Int {
        return Length(MissLocations(proposed, actual));
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
        let probabilities = EstimateClassificationProbabilities(
            tolerance, parameters, gates,
            Sampled(validationSchedule, features), nMeasurements
        );
        let localPL = InferredLabels(bias, probabilities);
        let nMismatches = NMismatches(localPL, labels);
        return ValidationResults(
            nMismatches
        );
    }

}
