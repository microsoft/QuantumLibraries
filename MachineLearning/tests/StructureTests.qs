// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.MachineLearning.Tests {
    open Microsoft.Quantum.Logical;
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Diagnostics;
    open Microsoft.Quantum.MachineLearning as ML;

    @Test("QuantumSimulator")
    function NQubitsRequiredFact() : Unit {
        let model = Default<ML.SequentialModel>()
            w/ Structure <- [
                ML.ControlledRotation((3, [7, 9]), PauliX, 0)
            ];
        let actual = ML.NQubitsRequired(model);
        EqualityFactI(actual, 10, "Wrong output from NQubitsRequired.");
    }

    function ExampleModel() : ML.SequentialModel {
        return Default<ML.SequentialModel>()
            w/ Structure <- [
                Default<ML.ControlledRotation>()
                    w/ TargetIndex <- 2
                    w/ ControlIndices <- [0]
                    w/ Axis <- PauliX
                    w/ ParameterIndex <- 0,
                Default<ML.ControlledRotation>()
                    w/ TargetIndex <- 0
                    w/ ControlIndices <- [1, 2]
                    w/ Axis <- PauliZ
                    w/ ParameterIndex <- 1
            ]
            w/ Parameters <- [
                1.234,
                2.345
            ];
    }

    operation ApplyExampleModelManually(register : Qubit[]) : Unit is Adj + Ctl {
        Controlled R([register[0]], (PauliX, 1.234, register[2]));
        Controlled R([register[1], register[2]], (PauliZ, 2.345, register[0]));
    }

    @Test("QuantumSimulator")
    operation TestApplySequentialClassifier() : Unit {
        AssertOperationsEqualReferenced(ML.NQubitsRequired(ExampleModel()),
            ML.ApplySequentialClassifier(ExampleModel(), _),
            ApplyExampleModelManually
        );
    }

    function EqualCR(x : ML.ControlledRotation, y : ML.ControlledRotation) : Bool {
        return x::Axis == y::Axis and
               All(EqualI, Zipped(x::ControlIndices, y::ControlIndices)) and
               x::TargetIndex == y::TargetIndex and
               x::ParameterIndex == y::ParameterIndex;
    }

    @Test("QuantumSimulator")
    function LocalRotationsLayerFact() : Unit {
        Fact(All(EqualCR, Zipped(
            ML.LocalRotationsLayer(3, PauliY),
            [
                Default<ML.ControlledRotation>()
                    w/ TargetIndex <- 0
                    w/ ControlIndices <- []
                    w/ Axis <- PauliY
                    w/ ParameterIndex <- 0,

                Default<ML.ControlledRotation>()
                    w/ TargetIndex <- 1
                    w/ ControlIndices <- []
                    w/ Axis <- PauliY
                    w/ ParameterIndex <- 1,

                Default<ML.ControlledRotation>()
                    w/ TargetIndex <- 2
                    w/ ControlIndices <- []
                    w/ Axis <- PauliY
                    w/ ParameterIndex <- 2
            ]
        )), "LocalRotationsLayer returned wrong output.");
    }

    @Test("QuantumSimulator")
    function PartialRotationsLayerFact() : Unit {
        Fact(All(EqualCR, Zipped(
            ML.PartialRotationsLayer([4, 5, 6], PauliY),
            [
                Default<ML.ControlledRotation>()
                    w/ TargetIndex <- 4
                    w/ ControlIndices <- []
                    w/ Axis <- PauliY
                    w/ ParameterIndex <- 0,

                Default<ML.ControlledRotation>()
                    w/ TargetIndex <- 5
                    w/ ControlIndices <- []
                    w/ Axis <- PauliY
                    w/ ParameterIndex <- 1,

                Default<ML.ControlledRotation>()
                    w/ TargetIndex <- 6
                    w/ ControlIndices <- []
                    w/ Axis <- PauliY
                    w/ ParameterIndex <- 2
            ]
        )), "PartialRotationsLayer returned wrong output.");
    }

    @Test("QuantumSimulator")
    function CyclicEntanglingLayerFact() : Unit {
        Fact(All(EqualCR, Zipped(
            ML.CyclicEntanglingLayer(3, PauliX, 2),
            [
                Default<ML.ControlledRotation>()
                    w/ TargetIndex <- 0
                    w/ ControlIndices <- [2]
                    w/ Axis <- PauliX
                    w/ ParameterIndex <- 0,

                Default<ML.ControlledRotation>()
                    w/ TargetIndex <- 1
                    w/ ControlIndices <- [0]
                    w/ Axis <- PauliX
                    w/ ParameterIndex <- 1,

                Default<ML.ControlledRotation>()
                    w/ TargetIndex <- 2
                    w/ ControlIndices <- [1]
                    w/ Axis <- PauliX
                    w/ ParameterIndex <- 2
            ]
        )), "CyclicEntanglingLayer returned wrong output.");
    }

    @Test("QuantumSimulator")
    function CombinedStructureFact() : Unit {
        let combined = ML.CombinedStructure([
            [
                Default<ML.ControlledRotation>()
                    w/ TargetIndex <- 0
                    w/ ControlIndices <- [2]
                    w/ Axis <- PauliX
                    w/ ParameterIndex <- 0,

                Default<ML.ControlledRotation>()
                    w/ TargetIndex <- 1
                    w/ ControlIndices <- [0]
                    w/ Axis <- PauliX
                    w/ ParameterIndex <- 1
            ],
            [
                Default<ML.ControlledRotation>()
                    w/ TargetIndex <- 2
                    w/ ControlIndices <- [1]
                    w/ Axis <- PauliZ
                    w/ ParameterIndex <- 0
            ]
        ]);
        Fact(All(EqualCR, Zipped(
            combined,
            [
                Default<ML.ControlledRotation>()
                    w/ TargetIndex <- 0
                    w/ ControlIndices <- [2]
                    w/ Axis <- PauliX
                    w/ ParameterIndex <- 0,

                Default<ML.ControlledRotation>()
                    w/ TargetIndex <- 1
                    w/ ControlIndices <- [0]
                    w/ Axis <- PauliX
                    w/ ParameterIndex <- 1,

                Default<ML.ControlledRotation>()
                    w/ TargetIndex <- 2
                    w/ ControlIndices <- [1]
                    w/ Axis <- PauliZ
                    w/ ParameterIndex <- 2
            ]
        )), "CombinedStructure returned wrong output.");
    }

}
