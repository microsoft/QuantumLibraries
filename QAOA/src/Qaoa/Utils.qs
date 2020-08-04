namespace Microsoft.Quantum.QAOA {

    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Measurement;


   operation MeasureAllAndReset(qubits: Qubit[]) : Bool[]
    {
        let N = Length(qubits);
        mutable results = new Bool[N];
        for (i in 0..N-1)
        {
            set results w/= i <- (MResetZ(qubits[i]) == One);
        }
        return results;
    }

    operation RunPhaseKickback(qubits: Qubit[], ancillaQubit: Qubit[], controlQubitsIndices: Int[], phaseExponent: Double) : Unit
    {
        for(i in 0..Length(controlQubitsIndices)-1)
        {
            CNOT(qubits[controlQubitsIndices[i]], ancillaQubit[0]);
		}

        R(PauliZ, phaseExponent, ancillaQubit[0]);

        for(i in 0..Length(controlQubitsIndices)-1)
        {
            CNOT(qubits[controlQubitsIndices[i]], ancillaQubit[0]);
		}
	}
}