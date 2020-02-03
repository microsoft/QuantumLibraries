// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.MachineLearning {
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Arithmetic;

    /// # Summary
    /// Describes a controlled rotation in terms of its target and control
    /// indices, rotation axis, and index into a model parameter vector.
    ///
    /// # Input
    /// ## TargetIndex
    /// Index of the target qubit for this controlled rotation.
    /// ## ControlIndices
    /// An array of the control qubit indices for this rotation.
    /// ## Axis
    /// The axis for this rotation.
    /// ## ParameterIndex
    /// An index into a model parameter vector describing the angle
    /// for this rotation.
    ///
    /// # Remarks
    /// An uncontrolled rotation can be represented by setting `ControlIndices`
    /// to an empty array of indexes, `new Int[0]`.
    ///
    /// # Example
    /// The following represents a rotation about the $X$-axis of the first
    /// qubit in a register, controlled on the second qubit, and with an
    /// angle given by the fourth parameter in a sequential model:
    /// ```Q#
    /// let controlledRotation = ControlledRotation(
    ///     (0, [1]),
    ///     PauliX,
    ///     3
    /// )
    /// ```
    newtype ControlledRotation = (
        (
            TargetIndex: Int,
            ControlIndices: Int[]
        ),
        Axis: Pauli,
        ParameterIndex: Int
    );

    /// # Summary
    /// Describes a quantum classifier model comprised of a sequence of
    /// parameterized and controlled rotations, an assignment of rotation
    /// angles, and a bias between the two classes recognized by the model.
    ///
    /// # Input
    /// ## Structure
    /// The sequence of controlled rotations used to classify inputs.
    /// ## Parameters
    /// An assignment of rotation angles to the given classification structure.
    /// ## Bias
    /// The bias between the two classes recognized by this classifier.
    ///
    /// # References
    /// - [arXiv:1804.00633](https://arxiv.org/abs/1804.00633)
    newtype SequentialModel = (
        Structure: ControlledRotation[],
        Parameters: Double[],
        Bias: Double
    );

    /// # Summary
    /// Describes an operation that prepares a given input to a sequential
    /// classifier.
    ///
    /// # Input
    /// ## NQubits
    /// The number of qubits on which the encoded input is defined.
    /// ## Prepare
    /// An operation which prepares the encoded input on a little-endian
    /// register of `NQubits` qubits.
    newtype StateGenerator = (
        NQubits: Int,
        Prepare: (LittleEndian => Unit is Adj + Ctl)
    );

    /// # Summary
    /// A sample, labeled with a class to which that sample belongs.
    ///
    /// # Input
    /// ## Features
    /// A vector of features for the given sample.
    /// ## Label
    /// An integer label for the class to which this sample belongs.
    newtype LabeledSample = (
        Features: Double[],
        Label: Int
    );

    /// # Summary
    /// A schedule for drawing batches from a set of samples.
    newtype SamplingSchedule = Range[];

    // Here, we define a couple private accessor functions for LabeledSample,
    // in lieu of having lambda support. These SHOULD NOT be used in external
    // code.
    function _Features(sample : LabeledSample) : Double[] { return sample::Features; }
    function _Label(sample : LabeledSample) : Int { return sample::Label; }

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
    /// The results from having validated a classifier against a set of
    /// samples.
    ///
    /// # Input
    /// ## NMisclassifications
    /// The number of misclassifications observed during validation.
    newtype ValidationResults = (
        NMisclassifications: Int
    );

    /// # Summary
    /// A collection of options to be used in training quantum classifiers.
    ///
    /// # Input
    /// ## LearningRate
    /// The learning rate by which gradients should be rescaled when updating
    /// model parameters during training steps.
    /// ## Tolerance
    /// The approximation tolerance to use when preparing samples as quantum
    /// states.
    /// ## MinibatchSize
    /// The number of samples to use in each training minibatch.
    /// ## NMeasurements
    /// The number of times to measure each classification result in order to
    /// estimate the classification probability.
    /// ## MaxEpochs
    /// The maximum number of epochs to train each model for.
    /// ## MaxStalls
    /// The maximum number of times a training epoch is allowed to stall
    /// (approximately zero gradient) before failing.
    /// ## StochasticRescaleFactor
    /// The amount to rescale stalled models by before retrying an update.
    ///
    /// # Remarks
    /// This UDT should not be created directly, but rather should be specified
    /// by calling @"microsoft.quantum.machinelearning.defaulttrainingoptions"
    /// and then using the `w/` operator to override different defaults.
    ///
    /// For example, to use 100,000 measurements:
    /// ```Q#
    /// let options = DefaultTrainingOptions() w/ NMeasurements <- 100000;
    /// ```
    ///
    /// # References
    /// - [arXiv:1804.00633](https://arxiv.org/abs/1804.00633)
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

}
