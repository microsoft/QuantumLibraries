// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {

	operation I(target : Qubit) : () {
		body {}
		adjoint self
		controlled auto
		controlled adjoint auto
	}

	operation NoOp(qubits: Qubit[]) : (){
		body {

		}
		adjoint auto
		controlled auto
		adjoint controlled auto
	}
	operation NoOp2(qubitsA: Qubit[], qubitsB: Qubit[]) : (){
		body {

		}
		adjoint auto
		controlled auto
		adjoint controlled auto
	}

	/// # Summary
	/// Ignores the output of an operation or function.
	function Ignore<'T>(value : 'T) : () {
		return ();
	}

}
