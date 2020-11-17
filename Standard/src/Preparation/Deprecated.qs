// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Preparation {
    open Microsoft.Quantum.Math;
    open Microsoft.Quantum.Arithmetic;

    /// # Summary
    /// Uses the Quantum ROM technique to represent a given density matrix.
    ///
    /// Given a list of $N$ coefficients $\alpha_j$, this returns a unitary $U$ that uses the Quantum-ROM
    /// technique to prepare
    /// an approximation  $\tilde\rho\sum_{j=0}^{N-1}p_j\ket{j}\bra{j}$ of the purification of the density matrix
    /// $\rho=\sum_{j=0}^{N-1}\frac{|alpha_j|}{\sum_k |\alpha_k|}\ket{j}\bra{j}$. In this approximation, the
    /// error $\epsilon$ is such that $|p_j-\frac{|alpha_j|}{\sum_k |\alpha_k|}|\le \epsilon / N$ and
    /// $\|\tilde\rho - \rho\| \le \epsilon$. In other words,
    /// $$
    /// \begin{align}
    /// U\ket{0}^{\lceil\log_2 N\rceil}\ket{0}^{m}=\sum_{j=0}^{N-1}\sqrt{p_j} \ket{j}\ket{\text{garbage}_j}.
    /// \end{align}
    /// $$
    ///
    /// # Input
    /// ## targetError
    /// The target error $\epsilon$.
    /// ## coefficients
    /// Array of $N$ coefficients specifying the probability of basis states.
    /// Negative numbers $-\alpha_j$ will be treated as positive $|\alpha_j|$.
    ///
    /// # Output
    /// ## First parameter
    /// A tuple `(x,(y,z))` where `x = y + z` is the total number of qubits allocated,
    /// `y` is the number of qubits for the `LittleEndian` register, and `z` is the Number
    /// of garbage qubits.
    /// ## Second parameter
    /// The one-norm $\sum_j |\alpha_j|$ of the coefficient array.
    /// ## Third parameter
    /// The unitary $U$.
    ///
    /// # Remarks
    /// ## Example
    /// The following code snippet prepares an purification of the $3$-qubit state
    /// $\rho=\sum_{j=0}^{4}\frac{|alpha_j|}{\sum_k |\alpha_k|}\ket{j}\bra{j}$, where
    /// $\vec\alpha=(1.0,2.0,3.0,4.0,5.0)$, and the error is `1e-3`;
    /// ```qsharp
    /// let coefficients = [1.0,2.0,3.0,4.0,5.0];
    /// let targetError = 1e-3;
    /// let ((nTotalQubits, (nIndexQubits, nGarbageQubits)), oneNorm, op) = QuantumROM(targetError, coefficients);
    /// using (indexRegister = Qubit[nIndexQubits]) {
    ///     using (garbageRegister = Qubit[nGarbageQubits]) {
    ///         op(LittleEndian(indexRegister), garbageRegister);
    ///     }
    /// }
    /// ```
    ///
    /// # References
    /// - Encoding Electronic Spectra in Quantum Circuits with Linear T Complexity
    ///   Ryan Babbush, Craig Gidney, Dominic W. Berry, Nathan Wiebe, Jarrod McClean, Alexandru Paler, Austin Fowler, Hartmut Neven
    ///   https://arxiv.org/abs/1805.03662
    @Deprecated("Microsoft.Quantum.Preparation.PurifiedMixedState")
    function QuantumROM(targetError: Double, coefficients: Double[])
    : ((Int, (Int, Int)), Double, ((LittleEndian, Qubit[]) => Unit is Adj + Ctl)) {
        let preparation = PurifiedMixedState(targetError, coefficients);
        return (
            preparation::Requirements!,
            preparation::Norm,
            IgnoreDataRegister(preparation::Prepare, _, _)
        );
    }

    internal operation IgnoreDataRegister(op : ((LittleEndian, Qubit[], Qubit[]) => Unit is Adj + Ctl), indexRegister : LittleEndian, garbageRegister : Qubit[]) : Unit is Adj + Ctl {
        op(indexRegister, new Qubit[0], garbageRegister);
    }

    /// # Summary
    /// Returns the total number of qubits that must be allocated
    /// to the operation returned by `QuantumROM`.
    ///
    /// # Input
    /// ## targetError
    /// The target error $\epsilon$.
    /// ## nCoeffs
    /// Number of coefficients specified in `QuantumROM`.
    ///
    /// # Output
    /// ## First parameter
    /// A tuple `(x,(y,z))` where `x = y + z` is the total number of qubits allocated,
    /// `y` is the number of qubits for the `LittleEndian` register, and `z` is the Number
    /// of garbage qubits.
    @Deprecated("Microsoft.Quantum.Preparation.PurifiedMixedStateRequirements")
    function QuantumROMQubitCount(targetError: Double, nCoeffs: Int)
    : (Int, (Int, Int)) {
        return (PurifiedMixedStateRequirements(targetError, nCoeffs))!;
    }

    /// # Summary
    /// Prepares a qubit in the +1 (`Zero`) eigenstate of the given Pauli operator.
    /// If the identity operator is given, then the qubit is prepared in the maximally
    /// mixed state.
    ///
    /// If the qubit was initially in the $\ket{0}$ state, this operation prepares the
    /// qubit in the $+1$ eigenstate of a given Pauli operator, or, for `PauliI`,
    /// in the maximally mixed state instead (see <xref:microsoft.quantum.preparation.preparesinglequbitidentity>).
    ///
    /// If the qubit was in a state other than $\ket{0}$, this operation applies the following gates:
    /// $H$ for `PauliX`, $HS$ for `PauliY`, $I$ for `PauliZ` and
    /// <xref:microsoft.quantum.preparation.preparesinglequbitidentity> for `PauliI`.
    ///
    /// # Input
    /// ## basis
    /// A Pauli operator $P$.
    /// ## qubit
    /// A qubit to be prepared.
    @Deprecated("Microsoft.Quantum.Preparation.PreparePauliEigenstate")
    operation PrepareQubit(basis : Pauli, qubit : Qubit) : Unit {
        PreparePauliEigenstate(basis, qubit);
    }

    /// # Summary
    /// Given a set of coefficients and a little-endian encoded quantum register,
    /// prepares an state on that register described by the given coefficients.
    ///
    /// # Description
    /// This operation prepares an arbitrary quantum
    /// state $\ket{\psi}$ with complex coefficients $r_j e^{i t_j}$ from
    /// the $n$-qubit computational basis state $\ket{0 \cdots 0}$.
    /// In particular, the action of this operation can be simulated by the
    /// a unitary transformation $U$ which acts on the all-zeros state as
    ///
    /// $$
    /// \begin{align}
    ///     U\ket{0...0}
    ///         & = \ket{\psi} \\\\
    ///         & = \frac{
    ///                 \sum_{j=0}^{2^n-1} r_j e^{i t_j} \ket{j}
    ///             }{
    ///                 \sqrt{\sum_{j=0}^{2^n-1} |r_j|^2}
    ///             }.
    /// \end{align}
    /// $$
    ///
    /// # Input
    /// ## coefficients
    /// Array of up to $2^n$ complex coefficients represented by their
    /// absolute value and phase $(r_j, t_j)$. The $j$th coefficient
    /// indexes the number state $\ket{j}$ encoded in little-endian format.
    ///
    /// ## qubits
    /// Qubit register encoding number states in little-endian format. This is
    /// expected to be initialized in the computational basis state
    /// $\ket{0...0}$.
    ///
    /// # Remarks
    /// Negative input coefficients $r_j < 0$ will be treated as though
    /// positive with value $|r_j|$. `coefficients` will be padded with
    /// elements $(r_j, t_j) = (0.0, 0.0)$ if fewer than $2^n$ are
    /// specified.
    ///
    /// # References
    /// - Synthesis of Quantum Logic Circuits
    ///   Vivek V. Shende, Stephen S. Bullock, Igor L. Markov
    ///   https://arxiv.org/abs/quant-ph/0406176
    ///
    /// # See Also
    /// - Microsoft.Quantum.Preparation.ApproximatelyPrepareArbitraryState
    @Deprecated("Microsoft.Quantum.Preparation.PrepareArbitraryStateCP")
    operation PrepareArbitraryState(coefficients : ComplexPolar[], qubits : LittleEndian) : Unit is Adj + Ctl {
        PrepareArbitraryStateCP(coefficients, qubits);
    }

    
    /// # Summary
    /// Returns an operation that prepares a specific quantum state.
    ///
    /// The returned operation $U$ prepares an arbitrary quantum
    /// state $\ket{\psi}$ with complex coefficients $r_j e^{i t_j}$ from
    /// the $n$-qubit computational basis state $\ket{0...0}$.
    ///
    /// The action of U on a newly-allocated register is given by
    /// $$
    /// \begin{align}
    ///     U\ket{0...0}=\ket{\psi}=\frac{\sum_{j=0}^{2^n-1}r_j e^{i t_j}\ket{j}}{\sqrt{\sum_{j=0}^{2^n-1}|r_j|^2}}.
    /// \end{align}
    /// $$
    ///
    /// # Input
    /// ## coefficients
    /// Array of up to $2^n$ complex coefficients represented by their
    /// absolute value and phase $(r_j, t_j)$. The $j$th coefficient
    /// indexes the number state $\ket{j}$ encoded in little-endian format.
    ///
    /// # Output
    /// A state-preparation unitary operation $U$.
    ///
    /// # Remarks
    /// Negative input coefficients $r_j < 0$ will be treated as though
    /// positive with value $|r_j|$. `coefficients` will be padded with
    /// elements $(r_j, t_j) = (0.0, 0.0)$ if fewer than $2^n$ are
    /// specified.
    ///
    /// ## Example
    /// The following snippet prepares the quantum state $\ket{\psi}=e^{i 0.1}\sqrt{1/8}\ket{0}+\sqrt{7/8}\ket{2}$
    /// in the qubit register `qubitsLE`.
    /// ```qsharp
    /// let amplitudes = [Sqrt(0.125), 0.0, Sqrt(0.875), 0.0];
    /// let phases = [0.1, 0.0, 0.0, 0.0];
    /// mutable complexNumbers = new ComplexPolar[4];
    /// for (idx in 0..3) {
    ///     set complexNumbers[idx] = ComplexPolar(amplitudes[idx], phases[idx]);
    /// }
    /// let op = StatePreparationComplexCoefficients(complexNumbers);
    /// using (qubits = Qubit[2]) {
    ///     let qubitsLE = LittleEndian(qubits);
    ///     op(qubitsLE);
    /// }
    /// ```
    @Deprecated("Microsoft.Quantum.Preparation.PrepareArbitraryStateCP")
    function StatePreparationComplexCoefficients(coefficients : ComplexPolar[]) : (LittleEndian => Unit is Adj + Ctl) {
        return PrepareArbitraryStateCP(coefficients, _);
    }

    
    /// # Summary
    /// Returns an operation that prepares the given quantum state.
    ///
    /// The returned operation $U$ prepares an arbitrary quantum
    /// state $\ket{\psi}$ with positive coefficients $\alpha_j\ge 0$ from
    /// the $n$-qubit computational basis state $\ket{0...0}$.
    ///
    /// The action of U on a newly-allocated register is given by
    /// $$
    /// \begin{align}
    ///     U \ket{0\cdots 0} = \ket{\psi} = \frac{\sum_{j=0}^{2^n-1}\alpha_j \ket{j}}{\sqrt{\sum_{j=0}^{2^n-1}|\alpha_j|^2}}.
    /// \end{align}
    /// $$
    ///
    /// # Input
    /// ## coefficients
    /// Array of up to $2^n$ coefficients $\alpha_j$. The $j$th coefficient
    /// indexes the number state $\ket{j}$ encoded in little-endian format.
    ///
    /// # Output
    /// A state-preparation unitary operation $U$.
    ///
    /// # Remarks
    /// Negative input coefficients $\alpha_j < 0$ will be treated as though
    /// positive with value $|\alpha_j|$. `coefficients` will be padded with
    /// elements $\alpha_j = 0.0$ if fewer than $2^n$ are specified.
    ///
    /// ## Example
    /// The following snippet prepares the quantum state $\ket{\psi}=\sqrt{1/8}\ket{0}+\sqrt{7/8}\ket{2}$
    /// in the qubit register `qubitsLE`.
    /// ```qsharp
    /// let amplitudes = [Sqrt(0.125), 0.0, Sqrt(0.875), 0.0];
    /// let op = StatePreparationPositiveCoefficients(amplitudes);
    /// using (qubits = Qubit[2]) {
    ///     let qubitsLE = LittleEndian(qubits);
    ///     op(qubitsLE);
    /// }
    /// ```
    @Deprecated("Microsoft.Quantum.Preparation.PrepareArbitraryStateD")
    function StatePreparationPositiveCoefficients(coefficients : Double[]) : (LittleEndian => Unit is Adj + Ctl) {
        return PrepareArbitraryStateD(coefficients, _);
    }


    /// # Summary
    /// Given a set of coefficients and a little-endian encoded quantum register,
    /// prepares an state on that register described by the given coefficients,
    /// up to a given approximation tolerance.
    ///
    /// # Description
    /// This operation prepares an arbitrary quantum
    /// state $\ket{\psi}$ with complex coefficients $r_j e^{i t_j}$ from
    /// the $n$-qubit computational basis state $\ket{0 \cdots 0}$.
    /// In particular, the action of this operation can be simulated by the
    /// a unitary transformation $U$ which acts on the all-zeros state as
    ///
    /// $$
    /// \begin{align}
    ///     U\ket{0...0}
    ///         & = \ket{\psi} \\\\
    ///         & = \frac{
    ///                 \sum_{j=0}^{2^n-1} r_j e^{i t_j} \ket{j}
    ///             }{
    ///                 \sqrt{\sum_{j=0}^{2^n-1} |r_j|^2}
    ///             }.
    /// \end{align}
    /// $$
    ///
    /// # Input
    /// ## tolerance
    /// The approximation tolerance to be used when preparing the given state.
    ///
    /// ## coefficients
    /// Array of up to $2^n$ complex coefficients represented by their
    /// absolute value and phase $(r_j, t_j)$. The $j$th coefficient
    /// indexes the number state $\ket{j}$ encoded in little-endian format.
    ///
    /// ## qubits
    /// Qubit register encoding number states in little-endian format. This is
    /// expected to be initialized in the computational basis state
    /// $\ket{0...0}$.
    ///
    /// # Remarks
    /// Negative input coefficients $r_j < 0$ will be treated as though
    /// positive with value $|r_j|$. `coefficients` will be padded with
    /// elements $(r_j, t_j) = (0.0, 0.0)$ if fewer than $2^n$ are
    /// specified.
    ///
    /// # References
    /// - Synthesis of Quantum Logic Circuits
    ///   Vivek V. Shende, Stephen S. Bullock, Igor L. Markov
    ///   https://arxiv.org/abs/quant-ph/0406176
    ///
    /// # See Also
    /// - Microsoft.Quantum.Preparation.ApproximatelyPrepareArbitraryState
    @Deprecated("Microsoft.Quantum.Preparation.ApproximatelyPrepareArbitraryStateCP")
    operation ApproximatelyPrepareArbitraryState(
        tolerance : Double,
        coefficients : ComplexPolar[],
        qubits : LittleEndian
    )
    : Unit is Adj + Ctl {
        ApproximatelyPrepareArbitraryStateCP(tolerance, coefficients, qubits);
    }
}
