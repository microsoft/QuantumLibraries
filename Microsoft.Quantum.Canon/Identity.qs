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

}