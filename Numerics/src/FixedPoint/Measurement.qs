// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Arithmetic {
    open Microsoft.Quantum.Convert;
    open Microsoft.Quantum.Measurement;
    open Microsoft.Quantum.Math;

    /// # Summary
    /// Measure a fixed-point number and return its value as Double.
    ///
    /// # Input
    /// ## fp
    /// Fixed-point number to measure.
    operation MeasureFixedPoint(fp : FixedPoint) : Double {
        let (p, xs) = fp!;
        let n = Length(xs);
        let sign = MResetZ(xs[n-1]) == One;
        mutable keepAdding = sign;
        mutable fpAsDouble = 0.;
        for (i in 0..n-2) {
            mutable currentRes = MResetZ(xs[i]) == (sign ? Zero | One);
            if (keepAdding) {
                set keepAdding = currentRes;
                set currentRes = not currentRes;
            }
            set fpAsDouble = fpAsDouble * 0.5 + (currentRes == true ? 1. | 0.);
        }
        return (sign ? -1.0 | 1.0) * fpAsDouble * PowD(2.0, IntAsDouble(p-2));
    }
}