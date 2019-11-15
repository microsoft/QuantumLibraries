// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.MachineLearning {
	open Microsoft.Quantum.Canon;
	open Microsoft.Quantum.Arithmetic;

	/// Qubit span of a multicontrolled single-qubit gate
	newtype GateSpan = (
		TargetIndex: Int,
		ControlIndices: Int[]
	);

	/// One-parameter controlled rotation gate triplet:
	/// (control structure, rotation axis, index of the rotation parameter)
	newtype ControlledRotation = (
        Span: GateSpan,
        Axis: Pauli,
        Index: Int
    );

	/// Abstraction for sequence of gates
	newtype GateSequence = ControlledRotation[];

	/// Abstraction for state preparation
	/// Fst(StateGenerator) is the number of qubits
	/// Snd(Stategenerator) is a circuit to prepare subject state
	newtype StateGenerator = (
        NQubits: Int,
        Apply: (LittleEndian => Unit is Adj + Ctl)
    );

	/// Convention: negative Snd(labledSample) signifies the last sample in a batch
	newtype LabeledSample = (
        Features: Double[],
        Label: Int
    );

	// Here, we define a couple private accessor functions for LabeledSample,
	// in lieu of having lambda support. These should not be used in external
	// code.
	function _Features(sample : LabeledSample) : Double[] { return sample::Features; }
	function _Label(sample : LabeledSample) : Int { return sample::Label; }

	/// Abstraction for a two-level range of indices
	newtype SamplingSchedule = Range[];

	newtype ValidationResults = (
		NMisclassifications: Int
	);



}
