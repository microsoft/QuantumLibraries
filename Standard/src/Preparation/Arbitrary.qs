// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Preparation {
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Arithmetic;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Math;
    open Microsoft.Quantum.Arrays;

    // This library returns operations that prepare a specified quantum state
    // from the computational basis state $\ket{0...0}$.

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
    function StatePreparationPositiveCoefficients (coefficients : Double[]) : (LittleEndian => Unit is Adj + Ctl) {
        let nCoefficients = Length(coefficients);
        mutable coefficientsComplexPolar = new ComplexPolar[nCoefficients];

        for (idx in 0 .. nCoefficients - 1) {
            set coefficientsComplexPolar w/= idx <- ComplexPolar(AbsD(coefficients[idx]), 0.0);
        }

        return PrepareArbitraryState(coefficientsComplexPolar, _);
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
    function StatePreparationComplexCoefficients (coefficients : ComplexPolar[]) : (LittleEndian => Unit is Adj + Ctl) {
        return PrepareArbitraryState(coefficients, _);
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
    operation PrepareArbitraryState(coefficients : ComplexPolar[], qubits : LittleEndian) : Unit is Adj + Ctl {
        ApproximatelyPrepareArbitraryState(0.0, coefficients, qubits);
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
    operation ApproximatelyPrepareArbitraryState(
        tolerance : Double,
        coefficients : ComplexPolar[],
        qubits : LittleEndian
    )
    : Unit is Adj + Ctl {
        (_CompileApproximateArbitraryStatePreparation(tolerance, coefficients, Length(qubits!)))(qubits);
    }

    operation _ApplyToLittleEndian(bareOp : ((Qubit[]) => Unit is Adj + Ctl), register : LittleEndian)
    : Unit is Adj + Ctl {
        bareOp(register!);
    }

    function _CompileApproximateArbitraryStatePreparation(
        tolerance : Double,
        coefficients : ComplexPolar[],
        nQubits : Int
    )
    : (LittleEndian => Unit is Adj + Ctl) {
        // pad coefficients at tail length to a power of 2.
        let coefficientsPadded = Padded(-2 ^ nQubits, ComplexPolar(0.0, 0.0), coefficients);
        let idxTarget = 0;
        let rngControl =
            // Determine what controls to apply to `op`.
            nQubits > 1
            ? (1 .. (nQubits - 1))
            | (1..0);
        let plan = _ApproximatelyUnprepareArbitraryStatePlan(
            tolerance, coefficientsPadded, (rngControl, idxTarget)
        );
        let unprepare = BoundCA(plan);
        return _ApplyToLittleEndian(Adjoint unprepare, _);
    }

    operation _ApplyMultiplexStep(
        tolerance : Double, disentangling : Double[], axis : Pauli,
        (rngControl : Range, idxTarget : Int),
        register : Qubit[]
    )
    : Unit is Adj + Ctl {
        let actualControl = LittleEndian(register[rngControl]);
        ApproximatelyMultiplexPauli(tolerance, disentangling, axis, actualControl, register[idxTarget]);
    }

    function _RangeLength(rng : Range) : Int {
        mutable len = 0;
        for (idx in rng) {
            set len += 1;
        }
        return len;
    }

    operation _ApplyGlobalRotationStep(
        angle : Double, idxTarget : Int, register : Qubit[]
    ) : Unit is Adj + Ctl {
        Exp([PauliI], angle, [register[idxTarget]]);
    }

    /// # Summary
    /// Implementation step of arbitrary state preparation procedure.
    ///
    /// # See Also
    /// - PrepareArbitraryState
    /// - Microsoft.Quantum.Canon.MultiplexPauli
    function _ApproximatelyUnprepareArbitraryStatePlan(
        tolerance : Double, coefficients : ComplexPolar[],
        (rngControl : Range, idxTarget : Int)
    )
    : (Qubit[] => Unit is Adj + Ctl)[] {
        mutable plan = new (Qubit[] => Unit is Adj + Ctl)[0];

        // For each 2D block, compute disentangling single-qubit rotation parameters
        let (disentanglingY, disentanglingZ, newCoefficients) = _StatePreparationSBMComputeCoefficients(coefficients);
        if (AnyOutsideToleranceD(tolerance, disentanglingZ)) {
            set plan += [_ApplyMultiplexStep(tolerance, disentanglingZ, PauliZ, (rngControl, idxTarget), _)];
        }
        if (AnyOutsideToleranceD(tolerance, disentanglingY)) {
            set plan += [_ApplyMultiplexStep(tolerance, disentanglingY, PauliY, (rngControl, idxTarget), _)];
        }

        // target is now in |0> state up to the phase given by arg of newCoefficients.

        // Continue recursion while there are control qubits.
        if (_RangeLength(rngControl) == 0) {
            let (abs, arg) = newCoefficients[0]!;
            if (AbsD(arg) > tolerance) {
                set plan += [_ApplyGlobalRotationStep(-1.0 * arg, idxTarget, _)];
            }
        } else {
            if (AnyOutsideToleranceCP(tolerance, newCoefficients)) {
                let newControl = (RangeStart(rngControl) + 1)..RangeStep(rngControl)..RangeEnd(rngControl);
                let newTarget = RangeStart(rngControl);
                set plan += _ApproximatelyUnprepareArbitraryStatePlan(tolerance, newCoefficients, (newControl, newTarget));
            }
        }

        return plan;
    }


    /// # Summary
    /// Computes the Bloch sphere coordinates for a single-qubit state.
    ///
    /// Given two complex numbers $a0, a1$ that represent the qubit state, computes coordinates
    /// on the Bloch sphere such that
    /// $a0 \ket{0} + a1 \ket{1} = r e^{it}(e^{-i \phi /2}\cos{(\theta/2)}\ket{0}+e^{i \phi /2}\sin{(\theta/2)}\ket{1})$.
    ///
    /// # Input
    /// ## a0
    /// Complex coefficient of state $\ket{0}$.
    /// ## a1
    /// Complex coefficient of state $\ket{1}$.
    ///
    /// # Output
    /// A tuple containing `(ComplexPolar(r, t), phi, theta)`.
    function BlochSphereCoordinates (a0 : ComplexPolar, a1 : ComplexPolar) : (ComplexPolar, Double, Double) {
        let abs0 = AbsComplexPolar(a0);
        let abs1 = AbsComplexPolar(a1);
        let arg0 = ArgComplexPolar(a0);
        let arg1 = ArgComplexPolar(a1);
        let r = Sqrt(abs0 * abs0 + abs1 * abs1);
        let t = 0.5 * (arg0 + arg1);
        let phi = arg1 - arg0;
        let theta = 2.0 * ArcTan2(abs1, abs0);
        return (ComplexPolar(r, t), phi, theta);
    }

    /// # Summary
    /// Implementation step of arbitrary state preparation procedure.
    /// # See Also
    /// - Microsoft.Quantum.Canon.PrepareArbitraryState
    function _StatePreparationSBMComputeCoefficients (coefficients : ComplexPolar[]) : (Double[], Double[], ComplexPolar[]) {
        mutable disentanglingZ = new Double[Length(coefficients) / 2];
        mutable disentanglingY = new Double[Length(coefficients) / 2];
        mutable newCoefficients = new ComplexPolar[Length(coefficients) / 2];

        for (idxCoeff in 0 .. 2 .. Length(coefficients) - 1) {
            let (rt, phi, theta) = BlochSphereCoordinates(coefficients[idxCoeff], coefficients[idxCoeff + 1]);
            set disentanglingZ w/= idxCoeff / 2 <- 0.5 * phi;
            set disentanglingY w/= idxCoeff / 2 <- 0.5 * theta;
            set newCoefficients w/= idxCoeff / 2 <- rt;
        }

        return (disentanglingY, disentanglingZ, newCoefficients);
    }

}


