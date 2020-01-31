namespace Microsoft.Quantum.MachineLearning {
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Logical;
    open Microsoft.Quantum.Canon;

    /// # Summary
    /// Given a set of inferred labels and a set of correct labels, returns
    /// indices for where each set of labels differs.
    ///
    /// # Input
    /// ## inferredLabels
    /// The labels inferred for a given training or validation set.
    /// ## actualLabels
    /// The true labels for a given training or validation set.
    ///
    /// # Output
    /// An array of indices `idx` such that
    /// `inferredLabels[idx] != actualLabels[idx]`.
    ///
    /// # Example
    /// ```Q#
    /// let misclassifications = Misclassifications([0, 1, 0, 0], [0, 1, 1, 0]);
    /// Message($"{misclassifications}"); // Will print [2].
    /// ```
    function Misclassifications(inferredLabels : Int[], actualLabels : Int[])
    : Int[] {
        return Where(
            NotEqualI,
            Zip(inferredLabels, actualLabels)
        );
    }

    /// # Summary
    /// Given a set of inferred labels and a set of correct labels, returns
    /// the number of indices at which each set of labels differ.
    ///
    /// # Input
    /// ## inferredLabels
    /// The labels inferred for a given training or validation set.
    /// ## actualLabels
    /// The true labels for a given training or validation set.
    ///
    /// # Output
    /// The number of indices `idx` such that
    /// `inferredLabels[idx] != actualLabels[idx]`.
    ///
    /// # Example
    /// ```Q#
    /// let nMisclassifications = NMisclassifications([1, 1, 0, 0], [0, 1, 1, 0]);
    /// Message($"{nMisclassifications}"); // Will print 2.
    /// ```
    function NMisclassifications(proposed: Int[], actual: Int[]): Int {
        return Length(Misclassifications(proposed, actual));
    }

    operation ValidateSequentialClassifier(
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
            tolerance, model,
            Sampled(validationSchedule, features), nMeasurements
        );
        let localPL = InferredLabels(model::Bias, probabilities);
        let nMisclassifications = NMisclassifications(localPL, labels);
        return ValidationResults(
            nMisclassifications
        );
    }

}
