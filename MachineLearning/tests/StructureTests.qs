// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.MachineLearning.Tests {
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
        EqualityFactI(actual, 10, "Wrong output from ScheduleLength.");
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

}
