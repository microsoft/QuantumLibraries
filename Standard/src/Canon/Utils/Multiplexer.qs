// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {
    open Microsoft.Quantum.Arithmetic;
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Convert;
    open Microsoft.Quantum.Diagnostics;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Math;

    /// # Summary
    /// Applies a Pauli rotation conditioned on an array of qubits.
    ///
    /// # Description
    /// This applies a multiply controlled unitary operation that performs
    /// rotations by angle $\theta_j$ about single-qubit Pauli operator $P$
    /// when controlled by the $n$-qubit number state $\ket{j}$.
    /// In particular, the action of this operation is represented by the
    /// unitary
    ///
    /// $$
    /// \begin{align}
    ///     U = \sum^{2^n - 1}_{j=0} \ket{j}\bra{j} \otimes e^{i P \theta_j}.
    /// \end{align}
    /// $$
    ///
    /// # Input
    /// ## coefficients
    /// Array of up to $2^n$ coefficients $\theta_j$. The $j$th coefficient
    /// indexes the number state $\ket{j}$ encoded in little-endian format.
    ///
    /// ## pauli
    /// Pauli operator $P$ that determines axis of rotation.
    ///
    /// ## control
    /// $n$-qubit control register that encodes number states $\ket{j}$ in
    /// little-endian format.
    ///
    /// ## target
    /// Single qubit register that is rotated by $e^{i P \theta_j}$.
    ///
    /// # Remarks
    /// `coefficients` will be padded with elements $\theta_j = 0.0$ if
    /// fewer than $2^n$ are specified.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApproximatelyMultiplexPauli
    operation MultiplexPauli(coefficients : Double[], pauli : Pauli, control : LittleEndian, target : Qubit)
    : Unit is Adj + Ctl {
        ApproximatelyMultiplexPauli(0.0, coefficients, pauli, control, target);
    }

    /// # Summary
    /// Applies a Pauli rotation conditioned on an array of qubits, truncating
    /// small rotation angles according to a given tolerance.
    ///
    /// # Description
    /// This applies a multiply controlled unitary operation that performs
    /// rotations by angle $\theta_j$ about single-qubit Pauli operator $P$
    /// when controlled by the $n$-qubit number state $\ket{j}$.
    /// In particular, the action of this operation is represented by the
    /// unitary
    ///
    /// $$
    /// \begin{align}
    ///     U = \sum^{2^n - 1}_{j=0} \ket{j}\bra{j} \otimes e^{i P \theta_j}.
    /// \end{align}
    /// ##
    ///
    /// # Input
    /// ## tolerance
    /// A tolerance below which small coefficients are truncated.
    ///
    /// ## coefficients
    /// Array of up to $2^n$ coefficients $\theta_j$. The $j$th coefficient
    /// indexes the number state $\ket{j}$ encoded in little-endian format.
    ///
    /// ## pauli
    /// Pauli operator $P$ that determines axis of rotation.
    ///
    /// ## control
    /// $n$-qubit control register that encodes number states $\ket{j}$ in
    /// little-endian format.
    ///
    /// ## target
    /// Single qubit register that is rotated by $e^{i P \theta_j}$.
    ///
    /// # Remarks
    /// `coefficients` will be padded with elements $\theta_j = 0.0$ if
    /// fewer than $2^n$ are specified.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.MultiplexPauli
    operation ApproximatelyMultiplexPauli(tolerance : Double, coefficients : Double[], pauli : Pauli, control : LittleEndian, target : Qubit)
    : Unit is Adj + Ctl {
        if pauli == PauliZ {
            let op = ApproximatelyMultiplexZ(tolerance, coefficients, control, _);
            op(target);
        } elif pauli == PauliX {
            let op = ApproximatelyMultiplexPauli(tolerance, coefficients, PauliZ, control, _);
            ApplyWithCA(H, op, target);
        } elif pauli == PauliY {
            let op = ApproximatelyMultiplexPauli(tolerance, coefficients, PauliX, control, _);
            ApplyWithCA(Adjoint S, op, target);
        } elif pauli == PauliI {
            ApproximatelyApplyDiagonalUnitary(tolerance, coefficients, control);
        } else {
            fail $"MultiplexPauli failed. Invalid pauli {pauli}.";
        }
    }

    /// # Summary
    /// Applies a Pauli Z rotation conditioned on an array of qubits.
    ///
    /// # Description
    /// This applies the multiply controlled unitary operation that performs
    /// rotations by angle $\theta_j$ about single-qubit Pauli operator $Z$
    /// when controlled by the $n$-qubit number state $\ket{j}$.
    /// In particular, this operation can be represented by the unitary
    ///
    /// $$
    /// \begin{align}
    ///     U = \sum^{2^n-1}_{j=0} \ket{j}\bra{j} \otimes e^{i Z \theta_j}.
    /// \end{align}
    /// $$
    ///
    /// # Input
    /// ## coefficients
    /// Array of up to $2^n$ coefficients $\theta_j$. The $j$th coefficient
    /// indexes the number state $\ket{j}$ encoded in little-endian format.
    ///
    /// ## control
    /// $n$-qubit control register that encodes number states $\ket{j}$ in
    /// little-endian format.
    ///
    /// ## target
    /// Single qubit register that is rotated by $e^{i P \theta_j}$.
    ///
    /// # Remarks
    /// `coefficients` will be padded with elements $\theta_j = 0.0$ if
    /// fewer than $2^n$ are specified.
    ///
    /// # References
    /// - Synthesis of Quantum Logic Circuits
    ///   Vivek V. Shende, Stephen S. Bullock, Igor L. Markov
    ///   https://arxiv.org/abs/quant-ph/0406176
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApproximatelyMultiplexZ
    operation MultiplexZ(coefficients : Double[], control : LittleEndian, target : Qubit)
    : Unit is Adj + Ctl {
        ApproximatelyMultiplexZ(0.0, coefficients, control, target);
    }

    internal function AnyOutsideToleranceD(tolerance : Double, coefficients : Double[]) : Bool {
        // NB: We don't currently use Any / Mapped for this, as we want to be
        //     able to short-circuit.
        for coefficient in coefficients {
            if AbsD(coefficient) >= tolerance {
                return true;
            }
        }
        return false;
    }

    internal function AnyOutsideToleranceCP(tolerance : Double, coefficients : ComplexPolar[]) : Bool {
        for coefficient in coefficients {
            if AbsComplexPolar(coefficient) > tolerance {
                return true;
            }
        }
        return false;
    }

    /// # Summary
    /// Applies a Pauli Z rotation conditioned on an array of qubits, truncating
    /// small rotation angles according to a given tolerance.
    ///
    /// # Description
    /// This applies the multiply controlled unitary operation that performs
    /// rotations by angle $\theta_j$ about single-qubit Pauli operator $Z$
    /// when controlled by the $n$-qubit number state $\ket{j}$.
    /// In particular, this operation can be represented by the unitary
    ///
    /// $$
    /// \begin{align}
    ///     U = \sum^{2^n-1}_{j=0} \ket{j}\bra{j} \otimes e^{i Z \theta_j}.
    /// \end{align}
    /// $$
    ///
    /// # Input
    /// ## tolerance
    /// A tolerance below which small coefficients are truncated.
    ///
    /// ## coefficients
    /// Array of up to $2^n$ coefficients $\theta_j$. The $j$th coefficient
    /// indexes the number state $\ket{j}$ encoded in little-endian format.
    ///
    /// ## control
    /// $n$-qubit control register that encodes number states $\ket{j}$ in
    /// little-endian format.
    ///
    /// ## target
    /// Single qubit register that is rotated by $e^{i P \theta_j}$.
    ///
    /// # Remarks
    /// `coefficients` will be padded with elements $\theta_j = 0.0$ if
    /// fewer than $2^n$ are specified.
    ///
    /// # References
    /// - Synthesis of Quantum Logic Circuits
    ///   Vivek V. Shende, Stephen S. Bullock, Igor L. Markov
    ///   https://arxiv.org/abs/quant-ph/0406176
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.MultiplexZ
    operation ApproximatelyMultiplexZ(tolerance : Double, coefficients : Double[], control : LittleEndian, target : Qubit) : Unit is Adj + Ctl {
        body (...) {
            // pad coefficients length at tail to a power of 2.
            let coefficientsPadded = Padded(-2 ^ Length(control!), 0.0, coefficients);

            if Length(coefficientsPadded) == 1 {
                // Termination case
                if AbsD(coefficientsPadded[0]) > tolerance {
                    Exp([PauliZ], coefficientsPadded[0], [target]);
                }
            } else {
                // Compute new coefficients.
                let (coefficients0, coefficients1) = MultiplexZCoefficients(coefficientsPadded);
                ApproximatelyMultiplexZ(tolerance, coefficients0, LittleEndian(Most(control!)), target);
                if AnyOutsideToleranceD(tolerance, coefficients1) {
                    within {
                        CNOT(Tail(control!), target);
                    } apply {
                        ApproximatelyMultiplexZ(tolerance,coefficients1, LittleEndian(Most(control!)), target);
                    }
                }
            }
        }

        controlled (controlRegister, ...) {
            // pad coefficients length to a power of 2.
            let coefficientsPadded = Padded(2 ^ (Length(control!) + 1), 0.0, Padded(-2 ^ Length(control!), 0.0, coefficients));
            let (coefficients0, coefficients1) = MultiplexZCoefficients(coefficientsPadded);
            ApproximatelyMultiplexZ(tolerance,coefficients0, control, target);
            if AnyOutsideToleranceD(tolerance,coefficients1) {
                within {
                    Controlled X(controlRegister, target);
                } apply {
                    ApproximatelyMultiplexZ(tolerance,coefficients1, control, target);
                }
            }
        }
    }

    /// # Summary
    /// Applies an array of complex phases to numeric basis states of a register
    /// of qubits.
    ///
    /// # Description
    /// This operation implements a diagonal unitary that applies a complex phase
    /// $e^{i \theta_j}$ on the $n$-qubit number state $\ket{j}$.
    /// In particular, this operation can be represented by the unitary
    ///
    /// $$
    /// \begin{align}
    ///     U = \sum^{2^n-1}_{j=0}e^{i\theta_j}\ket{j}\bra{j}.
    /// \end{align}
    /// $$
    ///
    /// # Input
    /// ## coefficients
    /// Array of up to $2^n$ coefficients $\theta_j$. The $j$th coefficient
    /// indexes the number state $\ket{j}$ encoded in little-endian format.
    /// ## qubits
    /// $n$-qubit control register that encodes number states $\ket{j}$ in
    /// little-endian format.
    ///
    /// # Remarks
    /// `coefficients` will be padded with elements $\theta_j = 0.0$ if
    /// fewer than $2^n$ are specified.
    ///
    /// # References
    /// - Synthesis of Quantum Logic Circuits
    ///   Vivek V. Shende, Stephen S. Bullock, Igor L. Markov
    ///   https://arxiv.org/abs/quant-ph/0406176
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApproximatelyApplyDiagonalUnitary
    operation ApplyDiagonalUnitary(coefficients : Double[], qubits : LittleEndian) : Unit is Adj + Ctl {
        ApproximatelyApplyDiagonalUnitary(0.0, coefficients, qubits);
    }

    /// # Summary
    /// Applies an array of complex phases to numeric basis states of a register
    /// of qubits, truncating small rotation angles according to a given
    /// tolerance.
    ///
    /// # Description
    /// This operation implements a diagonal unitary that applies a complex phase
    /// $e^{i \theta_j}$ on the $n$-qubit number state $\ket{j}$.
    /// In particular, this operation can be represented by the unitary
    ///
    /// $$
    /// \begin{align}
    ///     U = \sum^{2^n-1}_{j=0}e^{i\theta_j}\ket{j}\bra{j}.
    /// \end{align}
    /// $$
    ///
    /// # Input
    /// ## tolerance
    /// A tolerance below which small coefficients are truncated.
    ///
    /// ## coefficients
    /// Array of up to $2^n$ coefficients $\theta_j$. The $j$th coefficient
    /// indexes the number state $\ket{j}$ encoded in little-endian format.
    ///
    /// ## qubits
    /// $n$-qubit control register that encodes number states $\ket{j}$ in
    /// little-endian format.
    ///
    /// # Remarks
    /// `coefficients` will be padded with elements $\theta_j = 0.0$ if
    /// fewer than $2^n$ are specified.
    ///
    /// # References
    /// - Synthesis of Quantum Logic Circuits
    ///   Vivek V. Shende, Stephen S. Bullock, Igor L. Markov
    ///   https://arxiv.org/abs/quant-ph/0406176
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApplyDiagonalUnitary
    operation ApproximatelyApplyDiagonalUnitary(tolerance : Double, coefficients : Double[], qubits : LittleEndian)
    : Unit is Adj + Ctl {
        if IsEmpty(qubits!) {
            fail "operation ApplyDiagonalUnitary -- Number of qubits must be greater than 0.";
        }

        // pad coefficients length at tail to a power of 2.
        let coefficientsPadded = Padded(-2 ^ Length(qubits!), 0.0, coefficients);

        // Compute new coefficients.
        let (coefficients0, coefficients1) = MultiplexZCoefficients(coefficientsPadded);
        ApproximatelyMultiplexZ(tolerance,coefficients1, LittleEndian(Most(qubits!)), Tail(qubits!));

        if Length(coefficientsPadded) == 2 {
            // Termination case
            if AbsD(coefficients0[0]) > tolerance {
                Exp([PauliI], 1.0 * coefficients0[0], qubits!);
            }
        } else {
            ApproximatelyApplyDiagonalUnitary(tolerance, coefficients0, LittleEndian(Most(qubits!)));
        }
    }

    /// # Summary
    /// Implementation step of multiply-controlled Z rotations.
    /// # See Also
    /// - Microsoft.Quantum.Canon.MultiplexZ
    internal function MultiplexZCoefficients(coefficients : Double[]) : (Double[], Double[]) {
        let newCoefficientsLength = Length(coefficients) / 2;
        mutable coefficients0 = [0.0, size = newCoefficientsLength];
        mutable coefficients1 = [0.0, size = newCoefficientsLength];

        for idxCoeff in 0 .. newCoefficientsLength - 1 {
            set coefficients0 w/= idxCoeff <- 0.5 * (coefficients[idxCoeff] + coefficients[idxCoeff + newCoefficientsLength]);
            set coefficients1 w/= idxCoeff <- 0.5 * (coefficients[idxCoeff] - coefficients[idxCoeff + newCoefficientsLength]);
        }

        return (coefficients0, coefficients1);
    }

    /// # Summary
    /// Applies an array of operations controlled by an array of number states.
    ///
    /// That is, applies Multiply-controlled unitary operation $U$ that applies a
    /// unitary $V_j$ when controlled by $n$-qubit number state $\ket{j}$.
    ///
    /// $U = \sum^{2^n-1}_{j=0}\ket{j}\bra{j}\otimes V_j$.
    ///
    /// # Input
    /// ## unitaries
    /// Array of up to $2^n$ unitary operations. The $j$th operation
    /// is indexed by the number state $\ket{j}$ encoded in little-endian format.
    ///
    /// ## index
    /// $n$-qubit control register that encodes number states $\ket{j}$ in
    /// little-endian format.
    ///
    /// ## target
    /// Generic qubit register that $V_j$ acts on.
    ///
    /// # References
    /// - Encoding Electronic Spectra in Quantum Circuits with Linear T Complexity
    ///   Ryan Babbush, Craig Gidney, Dominic W. Berry, Nathan Wiebe, Jarrod McClean, Alexandru Paler, Austin Fowler, Hartmut Neven
    ///   https://arxiv.org/abs/1805.03662
    operation MultiplexOperations<'T>(unitaries : ('T => Unit is Adj + Ctl)[], index : LittleEndian, target : 'T)
    : Unit is Adj + Ctl {
        body (...) {
            let (N, n) = DimensionsForMultiplexer(Length(unitaries), index);

            if N == 1 { // base case
                Head(unitaries)(target);
            } else {
                let (most, tail) = MostAndTail(index![...n - 1]);
                let parts = Partitioned([2^(n - 1)], unitaries);

                within {
                    X(tail);
                } apply {
                    Controlled MultiplexOperations([tail], (parts[0], LittleEndian(most), target)); 
                }

                Controlled MultiplexOperations([tail], (parts[1], LittleEndian(most), target)); 
            }
        }

        controlled (ctls, ...) {
            let nCtls = Length(ctls);

            if nCtls == 0 {
                MultiplexOperations(unitaries, index, target);
            } elif nCtls == 1 {
                let (N, n) = DimensionsForMultiplexer(Length(unitaries), index);

                let ctl = Head(ctls);

                if N == 1 { // base case
                    Controlled (Head(unitaries))(ctls, target);
                } else {
                    use helper = Qubit();

                    let (most, tail) = MostAndTail(index![...n - 1]);
                    let parts = Partitioned([2^(n - 1)], unitaries);

                    within {
                        X(tail);
                    } apply {
                        ApplyAnd(ctl, tail, helper);
                    }

                    Controlled MultiplexOperations([helper], (parts[0], LittleEndian(most), target));

                    CNOT(ctl, helper);

                    Controlled MultiplexOperations([helper], (parts[1], LittleEndian(most), target));

                    Adjoint ApplyAnd(ctl, tail, helper);
                }
            } else {
                use helper = Qubit();
                within {
                    Controlled X(ctls, helper);
                } apply {
                    Controlled MultiplexOperations([helper], (unitaries, index, target));
                }
            }
        }
    }

    /// # Summary
    /// Validates and adjusts dimensions for address register
    ///
    /// # Description
    /// Given $N$ unitaries in `numUnitaries` and an address register of length $n'$,
    /// this function first checks whether $N \neq 0$ and $\lceil\log_2 N\rceil = n \le n'$,
    /// and then returns the tuple $(N, n)$.
    ///
    /// # Input
    /// ## numUnitaries
    /// The number of unitaries to multiplex.
    /// ## address
    /// The address register.
    internal function DimensionsForMultiplexer(numUnitaries : Int, address : LittleEndian) : (Int, Int) {
        let N = numUnitaries;
        Fact(N > 0, "data cannot be empty");

        let n = Ceiling(Lg(IntAsDouble(N)));
        Fact(Length(address!) >= n, $"address register is too small, requires at least {n} qubits");

        return (N, n);
    }

}
