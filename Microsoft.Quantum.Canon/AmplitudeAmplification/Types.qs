// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {
	open Microsoft.Quantum.Primitive;

	/// summary:
	///     Represents a reflection oracle.
	/// remarks:
	///     This oracle O = I - (1 - e^(- i ?)) |?><?| performs a partial
	///     reflection by a phase ? about a single pure state |?>.
	/// params:
	/// - The operation implementing the given reflection.
	///     - The phase ? by which to rotate the reflected subspace.
	///     - The qubit register on which to perform the given reflection>.
	newtype ReflectionOracle = ((Double, Qubit[]) => (): Adjoint, Controlled);

	///This oracle O|s>_a|?>_s = ? |t>_a U |?>_s + ... acts on the ancilla state |s>_a to implement the unitary U on any system state |?>_s with amplitude ? in the |t>_a basis.

	/// <summary>
	/// Represents an oracle for oblivious amplitude amplification.
	/// </summary>
	/// <remarks> 
	/// This oracle defined by O|s>_a|?>_s = ? |t>_a U |?>_s + ... acts on the ancilla state |s>_a to implement the unitary U on any system state |?>_s with amplitude ? in the |t>_a basis.
	/// The first parameter is the qubit register of |s>_a. The second parameter is the qubit register of |?>_s.
	/// </remarks>
	newtype ObliviousOracle = ((Qubit[], Qubit[]) => (): Adjoint, Controlled);

	/// <summary>
	/// Represents an oracle for state preparation.
	/// </summary>
	/// <remarks> 
	///This oracle defined by O|0>_f|0>_s = ? |1>_f |?>_s + ... acts on the on computational basis state |0>_f|0>_s to create the system state |?>_s with amplitude ? in the basis flagged by |1>_f.
	/// The first parameter is an index to the qubit register of |0>_f. The second parameter encompassed both registers.
	/// </remarks>
	newtype StateOracle = ((Int, Qubit[]) => (): Adjoint, Controlled);

	/// <summary>
	/// Represents an oracle for deterministic state preparation.
	/// </summary>
	/// <remarks> 
	///This oracle defined by O|0> = |?> acts on the on computational basis state |0> to create the state |?>.
	/// The first parameter is the qubit register of |?>.
	/// </remarks>
	newtype DeterministicStateOracle = (Qubit[] => (): Adjoint, Controlled);


	/// <summary>
	/// Phases for a sequence of partial reflections in amplitude amplification.
	/// </summary>
	/// <remarks> 
	/// The first parameter is an array of phases for reflection about the start state. The second parameter is an array of phases for reflection about the target state.
	/// Both arrays must be of equal length. Note that in many cases, the first phase about the start state and last phase about the target state introduces a global phase shift and may be set to 0. 
	/// </remarks>
	newtype AmpAmpReflectionPhases = (Double[], Double[]);

	/// <summary>
	/// Phases for a sequence of single-qubit rotations in amplitude amplification.
	/// </summary>
	/// <remarks> 
	/// The first parameter is an array of phases for reflections, expressed as a product of single-qubit rotations.
	/// [ G.H. Low, I. L. Chuang, https://arxiv.org/abs/1707.05391].
	/// </remarks>
	newtype AmpAmpRotationPhases = (Double[]);
	


}
