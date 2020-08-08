namespace Microsoft.Quantum.QAOA {

    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Measurement;

    /// # Summary
    /// Measures all qubits and resets them.
    ///
    /// # Input
    /// ## qubits
    /// Qubits to be measured.
    ///
    /// # Output
    /// Results of the measurement.

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

    /// # Summary
    /// Runs a phase kickback routine on an ancilla qubit given indices of control qubits and a phase.
    ///
    /// # Input
    /// ## qubits
    /// Qubits that are to aquire a phase.
    /// ## ancillaQubit
    /// Ancilla qubit.
    /// ## controlQubitsIndices
    /// List of indices of control qubits.
    /// ## phaseExponent
    /// Phase to be applied.

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