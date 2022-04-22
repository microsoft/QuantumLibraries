// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Arithmetic {
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Convert;
    open Microsoft.Quantum.Measurement;

    /// # Summary
    /// Measure a fixed-point number, returns its value as Double, and resets
    /// all the register to zero.
    ///
    /// # Input
    /// ## fp
    /// Fixed-point number to measure.
    operation MeasureFxP(fp : FixedPoint) : Double {
        let bits = ForEach(q => MResetZ(q) == One, fp::Register);
        return BoolArrayAsFixedPoint(fp::IntegerBits, bits);
    }
}
