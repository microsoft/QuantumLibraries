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
        let model = ML.SequentialModel([
            ML.ControlledRotation((3, [7, 9]), PauliX, 0),
            ML.ControlledRotation((8, []), PauliY, 1)
        ], [], 0.);
        let actual = ML.NQubitsRequired(model);
        EqualityFactI(actual, 10, "Wrong output from NQubitsRequired.");
    }

    function ExampleModel() : ML.SequentialModel {
        return ML.SequentialModel([
            ML.ControlledRotation((2, [0]), PauliX, 0),
            ML.ControlledRotation((0, [1, 2]), PauliZ, 1)],
            [1.234, 2.345],
            0.0);
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
                ML.ControlledRotation((0, []), PauliY, 0),
                ML.ControlledRotation((1, []), PauliY, 1),
                ML.ControlledRotation((2, []), PauliY, 2)
            ]
        )), "LocalRotationsLayer returned wrong output.");
    }

    @Test("QuantumSimulator")
    function PartialRotationsLayerFact() : Unit {
        Fact(All(EqualCR, Zipped(
            ML.PartialRotationsLayer([4, 5, 6], PauliY),
            [
                ML.ControlledRotation((4, []), PauliY, 0),
                ML.ControlledRotation((5, []), PauliY, 1),
                ML.ControlledRotation((6, []), PauliY, 2)
            ]
        )), "PartialRotationsLayer returned wrong output.");
    }

    @Test("QuantumSimulator")
    function CyclicEntanglingLayerFact() : Unit {
        Fact(All(EqualCR, Zipped(
            ML.CyclicEntanglingLayer(3, PauliX, 2),
            [
                ML.ControlledRotation((0, [2]), PauliX, 0),
                ML.ControlledRotation((1, [0]), PauliX, 1),
                ML.ControlledRotation((2, [1]), PauliX, 2)
            ]
        )), "CyclicEntanglingLayer returned wrong output.");
    }

    @Test("QuantumSimulator")
    function CombinedStructureFact() : Unit {
        let combined = ML.CombinedStructure([
            [
                ML.ControlledRotation((0, [2]), PauliX, 0),
                ML.ControlledRotation((1, [0]), PauliX, 1)
            ],
            [
                ML.ControlledRotation((2, [1]), PauliZ, 0)
            ]
        ]);
        Fact(All(EqualCR, Zipped(
            combined,
            [
                ML.ControlledRotation((0, [2]), PauliX, 0),
                ML.ControlledRotation((1, [0]), PauliX, 1),
                ML.ControlledRotation((2, [1]), PauliZ, 2)
            ]
        )), "CombinedStructure returned wrong output.");
    }

}
