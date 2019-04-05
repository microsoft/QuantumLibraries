namespace Microsoft.Quantum.Arithmetic {
    open Microsoft.Quantum.Primitive;
    open Microsoft.Quantum.Canon;

    /// # Summary
    /// This applies the in-place majority operation to 3 qubits.
    ///
    /// # Description
    /// If we denote the state of the target qubit as $\ket{z}$, and input states of
    /// the input qubits as $\ket{x}$ and $\ket{y}$, then
    /// this operation performs the following transformation:
    /// $\ket{xyz} \rightarrow \ket{x \oplus z} \ket{y \oplus z} \ket{\operatorname{MAJ} (x, y, z)}$.
    ///
    /// # Input
    /// ## input0
    /// The first input qubit.
    /// ## input1
    /// The second input qubit.
    /// ## target
    /// A qubit onto which the majority function will be applied.
    operation MAJ(input0 : Qubit, input1 : Qubit, target : Qubit) : Unit {
        body (...) {
            CNOT(target, input1);
            CNOT(target, input0);
            CCNOT(input1, input0, target);
        }
        adjoint auto;
        controlled auto;
        adjoint controlled auto;
    }

}
