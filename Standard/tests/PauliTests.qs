// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.
namespace Microsoft.Quantum.Tests {
    open Microsoft.Quantum.Preparation;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Measurement;
    open Microsoft.Quantum.Diagnostics;

    @Test("QuantumSimulator")
    operation TestMeasureWithScratch() : Unit {
        use register = Qubit[2];
        PrepareEntangledState([register[0]], [register[1]]);
        X(register[1]);
        let xxScratch = MeasureWithScratch([PauliX, PauliX], register);
        let xx = Measure([PauliX, PauliX], register);

        if xx != xxScratch {
            fail $"〈XX〉: MeasureWithScratch and Measure disagree";
        }

        let yyScratch = MeasureWithScratch([PauliY, PauliY], register);
        let yy = Measure([PauliY, PauliY], register);

        if yy != yyScratch {
            fail $"〈yy〉: MeasureWithScratch and Measure disagree";
        }

        let zzScratch = MeasureWithScratch([PauliZ, PauliZ], register);
        let zz = Measure([PauliZ, PauliZ], register);

        if zz != zzScratch {
            fail $"〈ZZ〉: MeasureWithScratch and Measure disagree";
        }

        ResetAll(register);
    }

}
