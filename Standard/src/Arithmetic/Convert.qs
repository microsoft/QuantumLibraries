// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Arithmetic {
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Arrays;

    /// # Summary
    /// Converts a `LittleEndian` qubit register to a `BigEndian` qubit
    /// register by reversing the qubit ordering.
    ///
    /// # Input
    /// ## input
    /// Qubit register in `LittleEndian` format.
    ///
    /// # Output
    /// Qubit register in `BigEndian` format.
    function LittleEndianAsBigEndian(input: LittleEndian) : BigEndian {
        return BigEndian(Reversed(input!));
    }

    /// # Summary
    /// Converts a `BigEndian` qubit register to a `LittleEndian` qubit
    /// register by reversing the qubit ordering.
    ///
    /// # Input
    /// ## input
    /// Qubit register in `BigEndian` format.
    ///
    /// # Output
    /// Qubit register in `LittleEndian` format.
    function BigEndianAsLittleEndian(input: BigEndian) : LittleEndian {
        return LittleEndian(Reversed(input!));
    }

}
