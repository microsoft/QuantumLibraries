// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Characterization {
    open Microsoft.Quantum.Measurement as Meas;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Arrays;

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Measurement.MeasureAllZ".
    @Deprecated("Microsoft.Quantum.Measurement.MeasureAllZ")
    operation MeasureAllZ(register : Qubit[]) : Result {
        return Measure(ConstantArray(Length(register), PauliZ), register);
    }

}
