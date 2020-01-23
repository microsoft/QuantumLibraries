// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

/// # Summary
/// This namespace contains functions and operations for preparing qubits into
/// different states, given registers in specific initial states.
///
/// # Description
/// This namespace provides operations for preparing different states, such as
/// entangled states, uniform superposition states, or arbitrary states
/// described by their coefficients in some basis.
///
/// By contrast with reset operations, the state preparation
/// operations provided by this namespace in general assume that quantum
/// registers provided as input are in fixed initial states, such as
/// the all-zeros state $\ket{00 \cdots 0}$, and prepare desired states
/// accordingly. Thus, preparation operations are more likely to be adjointable,
/// representing "unpreparation," useful for returning returning qubits
/// to their initial state without using intermediate measurements.
///
/// # See Also
/// - Microsoft.Quantum.Measurement
namespace Microsoft.Quantum.Preparation {}
