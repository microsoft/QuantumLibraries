// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.ErrorCorrection {
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Measurement;

    /// # Summary
    /// Syndrome measurement and the inverse of embedding.
    ///
    /// # Remarks
    /// $X$- and $Z$-stabilizers are not treated equally,
    /// which is due to the particular choice of the encoding circuit.
    /// This asymmetry leads to a different syndrome extraction routine.
    /// One could measure the syndrome by measuring multi-qubit Pauli operator
    /// directly on the code state, but for the distillation purpose
    /// the logical qubit is returned into a single qubit,
    /// in course of which the syndrome measurements can be done without further
    /// auxiliary qubits.
    ///
    /// Note that this operation is not marked as `internal`, as unit tests
    /// directly depend on this operation. As a future improvement, unit tests
    /// should be refactored to depend only on public callables directly.
    ///
    /// > [!WARNING]
    /// > This routine is tailored
    /// > to a particular encoding circuit for Steane's 7 qubit code;
    /// > if the encoding circuit is modified then the syndrome outcome
    /// > might have to be interpreted differently.
    ///
    /// # Output
    /// The logical qubit and a pair of integers for $X$-syndrome and $Z$-syndrome.
    /// They represent the index of the code qubit on which a single $X$- or $Z$-error
    /// would have caused the measured syndrome.
    operation _ExtractLogicalQubitFromSteaneCode(code : LogicalRegister)
    : (Qubit, Int, Int) {
        Adjoint SteaneCodeEncoderImpl((code!)[0 .. 0], (code!)[1 .. 6]);
        let x0 = M((code!)[6]);
        let x1 = M((code!)[1]);
        let x2 = M((code!)[3]);
        mutable xsyn = 0;

        if x0 == One {
            set xsyn = xsyn ^^^ 1;
        }

        if x1 == One {
            set xsyn = xsyn ^^^ 2;
        }

        if x2 == One {
            set xsyn = xsyn ^^^ 4;
        }

        set xsyn -= 1;

        // xsyn contains the qubit index (0..6) at which a single Z-error would
        // produce the given syndrome.
        let z0 = M((code!)[5]);
        let z1 = M((code!)[2]);
        let z2 = M((code!)[4]);
        mutable zsyn = 0;

        if z0 == One {
            set zsyn = zsyn ^^^ 1;
        }

        if z1 == One {
            set zsyn = zsyn ^^^ 2;
        }

        if z2 == One {
            set zsyn = zsyn ^^^ 5;
        }

        set zsyn -= 1;

        // zsyn contains the qubit index (0..6) at which a single X-error would
        // produce the given syndrome.
        return ((code!)[0], xsyn, zsyn);
    }


    /// # Summary
    /// Rotates a single qubit by π/4 about the Y-axis.
    ///
    /// # Description
    /// Performs a π/4 rotation about `Y`.
    ///
    /// The rotation is performed by consuming a magic
    /// state; that is, a copy of the state
    /// $$
    /// \begin{align}
    ///     \cos\frac{\pi}{8} \ket{0} + \sin \frac{\pi}{8} \ket{1}.
    /// \end{align}
    /// $$
    ///
    /// # Input
    /// ## data
    /// A qubit to be rotated about $Y$ by $\pi / 4$.
    ///
    /// ## magic
    /// A qubit initially in the magic state. Following application
    /// of this operation, `magic` is returned to the $\ket{0}$ state.
    ///
    /// # Remarks
    /// The following are equivalent:
    /// ```qsharp
    /// Ry(PI() / 4.0, data);
    /// ```
    /// and
    /// ```qsharp
    /// using (magic = Qubit()) {
    ///     Ry(PI() / 4.0, magic);
    ///     InjectPi4YRotation(data, magic);
    ///     Reset(magic);
    /// }
    /// ```
    ///
    /// This operation supports the `Adjoint` functor, in which
    /// case the same magic state is consumed, but the effect
    /// on the data qubit is a $-\pi/4$ $Y$-rotation.
    operation InjectPi4YRotation (data : Qubit, magic : Qubit)
    : Unit is Adj {
        body (...) {
            Adjoint S(data);
            CNOT(magic, data);
            S(data);
            let r = MResetY(magic);

            if r == One {
                // The following five operations are equivalent to
                // Ry( Pi()/2.0, data), up to global phase.
                // Since this operation does not support Controlled, we need
                // not worry about global phases.
                S(data);
                H(data);
                Adjoint S(data);
                H(data);
                Adjoint S(data);
            }
        }

        adjoint (...) {
            Adjoint S(data);
            CNOT(magic, data);
            S(data);
            let r = MResetY(magic);

            if r == Zero {
                S(data);
                H(data);
                S(data);
                H(data);
                Adjoint S(data);
            }
        }
    }


    /// # Summary
    /// Implements the Knill magic state distillation algorithm.
    ///
    /// # Description
    /// Given 15 approximate copies of a magic state
    /// $$
    /// \begin{align}
    ///     \cos\frac{\pi}{8} \ket{0} + \sin \frac{\pi}{8} \ket{1}
    /// \end{align},
    /// $$
    /// yields one higher-quality copy.
    ///
    /// # Input
    /// ## roughMagic
    /// A register of fifteen qubits containing approximate copies
    /// of a magic state. Following application of this distillation
    /// procedure, `roughMagic[0]` will contain one higher-quality
    /// copy, and the rest of the register will be reset to the
    /// $\ket{00\cdots 0}$ state.
    ///
    /// # Output
    /// If `true`, then the procedure succeeded and the higher-quality
    /// copy should be accepted. If `false`, the procedure failed, and
    /// the state of the register should be considered undefined.
    ///
    /// # Remarks
    /// We follow the algorithm of Knill.
    /// However, the present implementation is far from being optimal,
    /// as it uses too many qubits.
    /// The magic states are injected in this routine,
    /// in which case there are better protocols.
    ///
    /// # References
    /// - [Knill](https://arxiv.org/abs/quant-ph/0402171)
    operation KnillDistill (roughMagic : Qubit[]) : Bool {
        mutable accept = false;

        use scratch = Qubit[8];
        let aux = scratch[7];
        let code = scratch[0 .. 6];
        InjectPi4YRotation(code[0], roughMagic[14]);
        SteaneCodeEncoderImpl(code[0 .. 0], code[1 .. 6]);

        for idx in 0 .. 6 {
            Adjoint InjectPi4YRotation(code[idx], roughMagic[idx]);
            CNOT(code[idx], aux);
            InjectPi4YRotation(code[idx], roughMagic[idx + 7]);
        }

        let (logicalQubit, xsyn, zsyn) = _ExtractLogicalQubitFromSteaneCode(LogicalRegister(code));
        let m = M(aux);

        if (xsyn == -1 and zsyn == -1) and m == Zero {
            SWAP(logicalQubit, roughMagic[0]);
            set accept = true;
        }

        ResetAll(scratch);

    return accept;
    }

}
