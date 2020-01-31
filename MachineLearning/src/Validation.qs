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

    /// # Summary
    /// Validates a given sequential classifier against a given set of
    /// pre-labeled samples.
    ///
    /// # Input
    /// ## model
    /// The sequential model to be validated.
    /// ## samples
    /// The samples to be used to validate the given model.
    /// ## tolerance
    /// The approximation tolerance to use in encoding each sample as an input
    /// to the sequential classifier.
    /// ## nMeasurements
    /// The number of measurements to use in classifying each sample.
    /// ## validationSchedule
    /// The schedule by which samples should be drawn from the validation set.
    ///
    /// # Ouput
    /// The results of the given validation.
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
