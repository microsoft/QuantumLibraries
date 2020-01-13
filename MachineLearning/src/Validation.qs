namespace Microsoft.Quantum.MachineLearning {
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Logical;
    open Microsoft.Quantum.Canon;

    function Misclassifications(inferredLabels : Int[], actualLabels : Int[])
    : Int[] {
        return Where(
            NotEqualI,
            Zip(inferredLabels, actualLabels)
        );
    }

    function NMisclassifications(proposed: Int[], actual: Int[]): Int {
        return Length(Misclassifications(proposed, actual));
    }

    operation ValidateSequentialClassifier(
        gates: SequentialClassifierStructure,
        model : SequentialModel,
        samples : LabeledSample[],
        tolerance: Double,
        nMeasurements: Int,
        validationSchedule: SamplingSchedule
    )
    : ValidationResults {
        let features = Mapped(_Features, samples);
        let labels = Sampled(validationSchedule, Mapped(_Label, samples));
        let probabilities = EstimateClassificationProbabilities(
            tolerance, model::Parameters, gates,
            Sampled(validationSchedule, features), nMeasurements
        );
        let localPL = InferredLabels(model::Bias, probabilities);
        let nMisclassifications = NMisclassifications(localPL, labels);
        return ValidationResults(
            nMisclassifications
        );
    }

}
