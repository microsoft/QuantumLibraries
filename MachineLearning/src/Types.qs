// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.MachineLearning {
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Arithmetic;

    /// # Summary
    /// Represents the span of a single-qubit quantum operation with zero or
    /// more control qubits.
    newtype GateSpan = (
        TargetIndex: Int,
        ControlIndices: Int[]
    );

    /// # Summary
    /// Represents a controlled rotation of a given angle, about a given axis,
    /// and with a given target index and control register.
    newtype ControlledRotation = (
        Span: GateSpan,
        Axis: Pauli,
        Index: Int
    );

    /// Abstraction for sequence of gates
    newtype SequentialClassifierStructure = ControlledRotation[];

    /// # Summary
    /// Represents an operation used to prepare an input state on a register
    /// of qubits in the little-endian representation.
    ///
    /// # Input
    /// ## NQubits
    /// Specifies the number of qubits that this operation acts on.
    /// ## Apply
    /// Called to prepare the encoded input state on a given register.
    newtype StateGenerator = (
        NQubits: Int,
        Apply: (LittleEndian => Unit is Adj + Ctl)
    );

    /// # Summary
    /// Represents the features and classification label for a given sample.
    ///
    /// # Input
    /// ## Features
    /// The vector of features for the given sample.
    /// ## Label
    /// An integer label for the class to which the given sample belongs.
    newtype LabeledSample = (
        Features: Double[],
        Label: Int
    );

    // Here, we define a couple private accessor functions for LabeledSample,
    // in lieu of having lambda support. These should not be used in external
    // code.
    function _Features(sample : LabeledSample) : Double[] { return sample::Features; }
    function _Label(sample : LabeledSample) : Int { return sample::Label; }

    /// # Summary
    /// TODO
    /// Abstraction for a two-level range of indices
    newtype SamplingSchedule = Range[];

    /// # Summary
    /// Returns the number of elements in a given sampling schedule.
    ///
    /// # Input
    /// ## schedule
    /// A sampling schedule whose length is to be returned.
    ///
    /// # Output
    /// The number of elements in the given sampling schedule.
    function ScheduleLength(schedule : SamplingSchedule) : Int {
        mutable length = 0;
        for (range in schedule!) {
            for (index in range) {
                set length += 1;
            }
        }
        return length;
    }

    /// # Summary
    /// Samples a given array, using the given schedule.
    ///
    /// # Input
    /// ## schedule
    /// A schedule to use in sampling values.
    /// ## values
    /// An array of values to be sampled.
    ///
    /// # Output
    /// An array of elements from values, following the given schedule.
    function Sampled<'T>(schedule : SamplingSchedule, values : 'T[]) : 'T[] {
        mutable sampled = new 'T[0];
        for (range in schedule!) {
            for (index in range) {
                set sampled += [values[index]];
            }
        }
        return sampled;
    }

    /// # Summary
    /// Represents the results of validating a classifier against a validation
    /// set.
    ///
    /// # Input
    /// ## NMissclassifications
    /// The number of misclassifications observed during validation.
    newtype ValidationResults = (
        NMisclassifications: Int
    );

    /// # Summary
    /// Represents configurable options that can be set to control how a
    /// classifier is trained.
    ///
    /// # Input
    /// TODO
    newtype TrainingOptions = (
        LearningRate: Double,
        Tolerance: Double,
        MinibatchSize: Int,
        NMeasurements: Int,
        MaxEpochs: Int,
        MaxStalls: Int,
        StochasticRescaleFactor: Double
    );

    /// # Summary
    /// Returns a default set of options for training classifiers.
    ///
    /// # Output
    /// A reasonable set of default training options for use when training
    /// classifiers.
    ///
    /// # Example
    /// To use the default options, but with additional measurements, use the
    /// `w/` operator:
    /// ```Q#
    /// let options = DefaultTrainingOptions()
    ///     w/ NMeasurements <- 1000000;
    /// ```
    function DefaultTrainingOptions() : TrainingOptions {
        return TrainingOptions(
            0.1, 0.005, 15, 10000, 16, 8, 0.01
        );
    }

    newtype SequentialModel = (
        Parameters: Double[],
        Bias: Double
    );

}
