// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.MachineLearning {
    open Microsoft.Quantum.Math;
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Canon;

	/// What is the minimum number of qubits
	/// to support the subject gate sequence?
	/// Find the maximum qubit index m occuring
	/// in a gate sequence and return m+1
	function NQubitsRequired(seq : GateSequence) : Int {
		mutable nQubitsRequired = 0;
        for (gate in seq!) {
            set nQubitsRequired = Fold(
                MaxI, 0,
                gate::Span::ControlIndices + [
                    gate::Span::TargetIndex,
                    nQubitsRequired
                ]
            );
        }
        return nQubitsRequired;
	}

	/// Apply parameterized gate sequence to subject qubit register
	///
	operation _ApplyGates(parameters : Double[], gates: GateSequence, qubits : Qubit[]) : (Unit) is Adj + Ctl {
        //dumpRegisterToConsole(qubits);
        for (gate in gates!) {
            // let (gsp,p,ix) = gt!;
            if (gate::Index < Length(parameters)) {
                let input = (gate::Axis, parameters[gate::Index], qubits[gate::Span::TargetIndex]);
                if (IsEmpty(gate::Span::ControlIndices)) {
                    // Uncontrolled rotation of target
                    R(input);
                } else {
                    //TODO: should one validate the control indices first?
                    (Controlled R)(Subarray(gate::Span::ControlIndices, qubits), input);
                }
            }
        }
    }

	operation ApplyGates(parameters : Double[], gates: GateSequence): (Qubit[] => Unit is Adj + Ctl) {
        return _ApplyGates(parameters,gates,_);
	}

}
