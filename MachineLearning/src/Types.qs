// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.MachineLearning {
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Arithmetic;

    /// Qubit span of a multicontrolled single-qubit gate
    newtype GateSpan = (
        TargetIndex: Int,
        ControlIndices: Int[]
    );

    /// One-parameter controlled rotation gate triplet:
    /// (control structure, rotation axis, index of the rotation parameter)
    newtype ControlledRotation = (
        Span: GateSpan,
        Axis: Pauli,
        Index: Int
    );

    /// Abstraction for sequence of gates
    newtype SequentialClassifierStructure = ControlledRotation[];

    /// Abstraction for state preparation
    /// Fst(StateGenerator) is the number of qubits
    /// Snd(Stategenerator) is a circuit to prepare subject state
    newtype StateGenerator = (
        NQubits: Int,
        Apply: (LittleEndian => Unit is Adj + Ctl)
    );

    /// Convention: negative Snd(labledSample) signifies the last sample in a batch
    newtype LabeledSample = (
        Features: Double[],
        Label: Int
    );

    // Here, we define a couple private accessor functions for LabeledSample,
    // in lieu of having lambda support. These should not be used in external
    // code.
    function _Features(sample : LabeledSample) : Double[] { return sample::Features; }
    function _Label(sample : LabeledSample) : Int { return sample::Label; }

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

    newtype ValidationResults = (
        NMisclassifications: Int
    );

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
