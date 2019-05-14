// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Arithmetic {
    open Microsoft.Quantum.Convert;
    open Microsoft.Quantum.Primitive;
    open Microsoft.Quantum.Math;

    /// # Summary
    /// Initialize a quantum fixed-point number to a classical constant.
    ///
    /// # Input
    /// ## constant
    /// Constant to which to initialize the quantum fixed-point number.
    /// ## fp
    /// Fixed-point number (of type FixedPoint) to initialize.
    operation FixedPointInit(constant : Double, fp : FixedPoint) : Unit{
        body (...) {
            let (p, q) = fp!;
            let n = Length(q);
            let sign = constant < 0.;
            mutable rescaledConstant = PowD(2., IntAsDouble(n-p)) * AbsD(constant) + 0.5;
            mutable keepAdding = sign;
            for (i in 0..n-1) {
                let intConstant = Floor(rescaledConstant);
                set rescaledConstant = 0.5 * rescaledConstant;
                mutable currentBit = (intConstant &&& 1) == (sign ? 0 | 1);
                if (keepAdding) {
                    set keepAdding = currentBit;
                    set currentBit = not currentBit;
                }
                if (currentBit) {
                    X(q[i]);
                }
            }
        }
        controlled auto;
        adjoint self;
        adjoint controlled auto;
    }
}