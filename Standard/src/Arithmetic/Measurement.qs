// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Arithmetic {
    open Microsoft.Quantum.Measurement;
    open Microsoft.Quantum.Convert;

    /// # Summary
    /// Measures the content of a quantum register and converts
    /// it to an integer. The measurement is performed with respect
    /// to the standard computational basis, i.e., the eigenbasis of `PauliZ`.
    ///
    /// # Input
    /// ## target
    /// A quantum register in the little-endian encoding.
    ///
    /// # Output
    /// An unsigned integer that contains the measured value of `target`.
    ///
    /// # Remarks
    /// This operation resets its input register to the $\ket{00\cdots 0}$ state,
    /// suitable for releasing back to a target machine.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.MeasureIntegerBE
    operation MeasureInteger(target : LittleEndian) : Int {
        mutable results = new Result[Length(target!)];

        for (idx in 0 .. Length(target!) - 1) {
            set results[idx] = MResetZ((target!)[idx]);
        }

        return ResultArrayAsInt(results);
    }

}
