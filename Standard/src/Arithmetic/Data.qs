// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Arithmetic {

    /// # Summary
    /// Register that encodes an unsigned integer in little-endian order. The
    /// qubit with index `0` encodes the lowest bit of an unsigned integer.
    ///
    /// # Remarks
    /// We abbreviate `LittleEndian` as `LE` in the documentation.
    newtype LittleEndian = Qubit[];

    /// # Summary
    /// Register that encodes an unsigned integer in big-endian order. The
    /// qubit with index `0` encodes the highest bit of an unsigned integer.
    ///
    /// # Remarks
    /// We abbreviate `BigEndian` as `BE` in the documentation.
    newtype BigEndian = Qubit[];

    /// # Summary
    /// Little-endian unsigned integers in QFT basis.
	///
    /// For example, if $\ket{x}$ is the little-endian encoding of the integer
    /// $x$ in the computational basis,
    /// then $\operatorname{QFTLE} \ket{x}$ is the encoding of $x$ in the QFT
    /// basis.
    ///
    /// # Remarks
    /// We abbreviate `PhaseLittleEndian` as `PhaseLE` in the documentation.
    ///
    /// # See Also
    /// - QFT
    /// - QFTLE
    newtype PhaseLittleEndian = Qubit[];

}
