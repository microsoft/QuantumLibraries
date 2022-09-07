// Copyright (c) Microsoft Corporation.
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
    operation PrepareArbitraryStateCP(coefficients : ComplexPolar[], qubits : LittleEndian) : Unit is Adj + Ctl {
        ApproximatelyPrepareArbitraryStateCP(0.0, coefficients, qubits);
    }

    /// # Summary
    /// Given a set of coefficients and a little-endian encoded quantum register,
    /// prepares an state on that register described by the given coefficients.
    ///
    /// # Description
    /// This operation prepares an arbitrary quantum
    /// state $\ket{\psi}$ with positive coefficients $\alpha_j\ge 0$ from
    /// the $n$-qubit computational basis state $\ket{0...0}$.
    ///
    /// The action of U on the all-zeros state is given by
    /// $$
    /// \begin{align}
    ///     U \ket{0\cdots 0} = \ket{\psi} = \frac{\sum_{j=0}^{2^n-1}\alpha_j \ket{j}}{\sqrt{\sum_{j=0}^{2^n-1}|\alpha_j|^2}}.
    /// \end{align}
    /// $$
    ///
    /// # Input
    /// ## coefficients
    /// Array of up to $2^n$ real coefficients. The $j$th coefficient
    /// indexes the number state $\ket{j}$ encoded in little-endian format.
    ///
    /// ## qubits
    /// Qubit register encoding number states in little-endian format. This is
    /// expected to be initialized in the computational basis state
    /// $\ket{0...0}$.
    ///
    /// # Remarks
    /// Negative input coefficients $\alpha_j < 0$ will be treated as though
    /// positive with value $|\alpha_j|$. `coefficients` will be padded with
    /// elements $\alpha_j = 0.0$ if fewer than $2^n$ are specified.
    ///
    /// # Example
    /// The following snippet prepares the quantum state $\ket{\psi}=\sqrt{1/8}\ket{0}+\sqrt{7/8}\ket{2}$
    /// in the qubit register `qubitsLE`.
    /// ```qsharp
    /// let amplitudes = [Sqrt(0.125), 0.0, Sqrt(0.875), 0.0];
    /// use qubits = Qubit[2];
    /// let qubitsLE = LittleEndian(qubits);
    /// PrepareArbitraryStateD(amplitudes, qubitsLE);
    /// ```
    ///
    /// # References
    /// - Synthesis of Quantum Logic Circuits
    ///   Vivek V. Shende, Stephen S. Bullock, Igor L. Markov
    ///   https://arxiv.org/abs/quant-ph/0406176
    ///
    /// # See Also
    /// - Microsoft.Quantum.Preparation.ApproximatelyPrepareArbitraryState
    operation PrepareArbitraryStateD(coefficients : Double[], qubits : LittleEndian) : Unit is Adj + Ctl {
        ApproximatelyPrepareArbitraryStateD(0.0, coefficients, qubits);
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
    operation ApproximatelyPrepareArbitraryStateCP(
        tolerance : Double,
        coefficients : ComplexPolar[],
        qubits : LittleEndian
    )
    : Unit is Adj + Ctl {
        (_CompileApproximateArbitraryStatePreparation(tolerance, coefficients, Length(qubits!)))(qubits);
    }

    /// # Summary
    /// Given a set of coefficients and a little-endian encoded quantum register
    /// of unentangled qubits, all of which are in zero state, prepares a state
    /// on that register described by the given coefficients. Uses state emulation
    /// if supported by the target.
    ///
    /// # Notes
    /// If the register isn't in zero state, the emulation will fail and fallback
    /// to quantum state preparation.
    /// This operation doesn't provide Adj/Ctr variants, because, in general, there
    /// are no efficient emulation algorithms for those.
    ///
    /// For internal use only, until proposal https://github.com/microsoft/qsharp-language/pull/41
    /// is finalized and implemented.
    operation _PrepareAmplitudesFromZeroState(coefficients : ComplexPolar[], qubits : LittleEndian) : Unit {
        ApproximatelyPrepareArbitraryStateCP(0.0, coefficients, qubits);
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
    /// Array of up to $2^n$ real coefficients. The $j$th coefficient
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
    operation ApproximatelyPrepareArbitraryStateD(
        tolerance : Double,
        coefficients : Double[],
        qubits : LittleEndian
    )
    : Unit is Adj + Ctl {
        let coefficientsAsComplexPolar = Mapped(Compose(ComplexPolar(_, 0.0), AbsD), coefficients);
        ApproximatelyPrepareArbitraryStateCP(tolerance, coefficientsAsComplexPolar, qubits);
    }

    /// # Summary
    /// Applies an operation to the underlying qubits making up a little-endian
    /// register. This operation is marked as internal, as a little-endian
    /// register is intended to be "opaque," such that this is an implementation
    /// detail only.
    internal operation ApplyToLittleEndian(bareOp : ((Qubit[]) => Unit is Adj + Ctl), register : LittleEndian)
    : Unit is Adj + Ctl {
        bareOp(register!);
    }

    // NB: This is currently not marked as internal, as the QML library
    //     currently uses this function. Please see the relevant GitHub issue
    //     at https://github.com/microsoft/QuantumLibraries/issues/239.
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
        let plan = ApproximatelyUnprepareArbitraryStatePlan(
            tolerance, coefficientsPadded, (rngControl, idxTarget)
        );
        let unprepare = BoundCA(plan);
        return ApplyToLittleEndian(Adjoint unprepare, _);
    }

    internal operation ApplyMultiplexStep(
        tolerance : Double, disentangling : Double[], axis : Pauli,
        (rngControl : Range, idxTarget : Int),
        register : Qubit[]
    )
    : Unit is Adj + Ctl {
        let actualControl = LittleEndian(register[rngControl]);
        ApproximatelyMultiplexPauli(tolerance, disentangling, axis, actualControl, register[idxTarget]);
    }

    internal operation ApplyGlobalRotationStep(
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
    internal function ApproximatelyUnprepareArbitraryStatePlan(
        tolerance : Double, coefficients : ComplexPolar[],
        (rngControl : Range, idxTarget : Int)
    )
    : (Qubit[] => Unit is Adj + Ctl)[] {
        mutable plan = [];

        // For each 2D block, compute disentangling single-qubit rotation parameters
        let (disentanglingY, disentanglingZ, newCoefficients) = StatePreparationSBMComputeCoefficients(coefficients);
        if (AnyOutsideToleranceD(tolerance, disentanglingZ)) {
            set plan += [ApplyMultiplexStep(tolerance, disentanglingZ, PauliZ, (rngControl, idxTarget), _)];
        }
        if (AnyOutsideToleranceD(tolerance, disentanglingY)) {
            set plan += [ApplyMultiplexStep(tolerance, disentanglingY, PauliY, (rngControl, idxTarget), _)];
        }

        // target is now in |0> state up to the phase given by arg of newCoefficients.

        // Continue recursion while there are control qubits.
        if (IsRangeEmpty(rngControl)) {
            let (abs, arg) = newCoefficients[0]!;
            if (AbsD(arg) > tolerance) {
                set plan += [ApplyGlobalRotationStep(-1.0 * arg, idxTarget, _)];
            }
        } elif (AnyOutsideToleranceCP(tolerance, newCoefficients)) {
            let newControl = (RangeStart(rngControl) + 1)..RangeStep(rngControl)..RangeEnd(rngControl);
            let newTarget = RangeStart(rngControl);
            set plan += ApproximatelyUnprepareArbitraryStatePlan(tolerance, newCoefficients, (newControl, newTarget));
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
    /// - Microsoft.Quantum.Preparation.PrepareArbitraryState
    internal function StatePreparationSBMComputeCoefficients(coefficients : ComplexPolar[]) : (Double[], Double[], ComplexPolar[]) {
        mutable disentanglingZ = [0.0, size = Length(coefficients) / 2];
        mutable disentanglingY = [0.0, size = Length(coefficients) / 2];
        mutable newCoefficients = [ComplexPolar(0.0, 0.0), size = Length(coefficients) / 2];

        for idxCoeff in 0 .. 2 .. Length(coefficients) - 1 {
            let (rt, phi, theta) = BlochSphereCoordinates(coefficients[idxCoeff], coefficients[idxCoeff + 1]);
            set disentanglingZ w/= idxCoeff / 2 <- 0.5 * phi;
            set disentanglingY w/= idxCoeff / 2 <- 0.5 * theta;
            set newCoefficients w/= idxCoeff / 2 <- rt;
        }

        return (disentanglingY, disentanglingZ, newCoefficients);
    }

}
