// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Arithmetic {
    open Microsoft.Quantum.Primitive;
    open Microsoft.Quantum.Canon;

    /// # Summary
    /// Converts a `LittleEndian` qubit register to a `BigEndian` Qubit
    /// register by reversing the qubit ordering.
    ///
    /// # Input
    /// ## input
    /// Qubit register in `LittleEndian` format.
    ///
    /// # Output
    /// Qubit register in `BigEndian` format.
    function LittleEndianToBigEndian(input: LittleEndian) : BigEndian {
        return BigEndian(Reverse(input!));
    }

    /// # Summary
    /// Converts a `BigEndian` qubit register to a `LittleEndian` Qubit
    /// register by reversing the qubit ordering.
    ///
    /// # Input
    /// ## input
    /// Qubit register in `BigEndian` format.
    ///
    /// # Output
    /// Qubit register in `LittleEndian` format.
    function BigEndianToLittleEndian(input: BigEndian) : LittleEndian {
        return LittleEndian(Reverse(input!));
    }

}
