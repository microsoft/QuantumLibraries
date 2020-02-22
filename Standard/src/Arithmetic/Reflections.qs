namespace Microsoft.Quantum.Arithmetic {
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Convert;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Logical;

	/// # Summary
    /// Reflects a quantum register about a given classical integer.
    ///
    /// # Description
    /// Given a quantum register initially in the state $\sum_i \alpha_i \ket{i}$,
    /// where each $\ket{i}$ is a basis state representing an integer $i$,
    /// reflects the state of the register about the basis state for a given
    /// integer $\ket{j}$,
    /// $$
    ///     \sum_i (-1)^{ \delta_{ij} } \alpha_i \ket{i}
    /// $$
    ///
    /// # Input
    /// ## index
    /// The classical integer indexing the basis state about which to reflect.
    ///
    /// # Remarks
    /// This operation is implemented in-place, without explicit allocation of
    /// additional auxillary qubits.
	operation ReflectAboutInteger(index : Int, reg : LittleEndian) : Unit is Adj + Ctl {
        within {
            // We want to reduce to the problem of reflecting about the all-ones
            // state. To do that, we apply our reflection within an application
            // of X instructions that flip all the zeros in our index.
            ApplyToEachCA(
                CControlledCA(X),
                Zip(Mapped(Not, IntAsBoolArray(index, Length(reg!))), reg!)
            );
        } apply {
            Controlled Z(Most(reg!), Tail(reg!));
        }
	}

}
