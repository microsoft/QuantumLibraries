// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Characterization {
    open Microsoft.Quantum.Measurement as Meas;

    /// # Deprecated
    /// Please use @"microsoft.quantum.measurement.measureAllZ".
    @Deprecated("Microsoft.Quantum.Measurement.MeasureAllZ")
    operation MeasureAllZ(register : Qubit[]) : Result {
        return Meas.MeasureAllZ(register);
    }

}
