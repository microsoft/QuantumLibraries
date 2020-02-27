// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.MachineLearning.Tests {
    open Microsoft.Quantum.Diagnostics;
    open Microsoft.Quantum.MachineLearning as ML;

    @Test("QuantumSimulator")
    function ScheduleLengthFact() : Unit {
        let actualLength = ML.ScheduleLength(ML.SamplingSchedule([0..4, 1..2..5]));
        EqualityFactI(actualLength, 5 + 3, "Wrong output from ScheduleLength.");
    }

    @Test("QuantumSimulator")
    function SampledFact() : Unit {
        let actuallySampled = ML.Sampled(ML.SamplingSchedule([0..4, 1..2..5]), [0, 10, 20, 30, 40, 50, 60, 70]);
        AllEqualityFactI(actuallySampled, [0, 10, 20, 30, 40, 10, 30, 50], "Wrong output from Sampled.");
    }

}
