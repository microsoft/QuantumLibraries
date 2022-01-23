// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Chemistry.JordanWigner {
    open Microsoft.Quantum.Arithmetic;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Chemistry;
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Math;


    /// # Summary
    /// Applies a sequence of Z operations and either an X or Y operation to
    /// a register of qubits, where the selection of target qubits and basis
    /// are conditioned on the state of a control register.
    ///
    /// # Description
    /// This operation can be described by a unitary matrix $U$ that applies
    /// the Pauli string on $(X^{z+1}\_pY^{z}\_p)Z\_{p-1}...Z_0$ on
    /// qubits $0..p$ conditioned on an index $z\in\{0,1\}$ and $p$.
    ///
    /// That is,
    /// $$
    /// \begin{align}
    /// U\ket{z}\ket{p}\ket{\psi} = \ket{z}\ket{p}(X^{z+1}\_pY^{z}\_p)Z\_{p-1}...Z_0\ket{\psi}
    /// \end{align}
    /// $$
    ///
    /// # Input
    /// ## pauliBasis
    /// When this qubit is in state $\ket{0}$, an `X` operation is applied. When it is in state $\ket{1}$, `Y` is applied.
    /// ## indexRegister
    /// The state $\ket{p}$ of this register determines the qubit on which `X` or `Y` is applied.
    /// ## targetRegister
    /// Register of qubits on which the Pauli operators are applied.
    ///
    /// # References
    /// - Encoding Electronic Spectra in Quantum Circuits with Linear T Complexity
    ///   Ryan Babbush, Craig Gidney, Dominic W. Berry, Nathan Wiebe, Jarrod McClean, Alexandru Paler, Austin Fowler, Hartmut Neven
    ///   https://arxiv.org/abs/1805.03662
    operation OptimizedBEXY (pauliBasis : Qubit, indexRegister : LittleEndian, targetRegister : Qubit[])
    : Unit is Adj + Ctl {
        let unitaryGenerator = (Length(targetRegister), _OptimizedBEXY_);

        use accumulator = Qubit();

        // This assumes that MultiplexOperationsFromGenerator applies unitaries indexed in unitaryGenerator in ascending order.
        X(accumulator);
        MultiplexOperationsFromGenerator(unitaryGenerator, indexRegister, (pauliBasis, accumulator, targetRegister));
        // If indexRegister encodes an integer that is larger than Length(targetRegister),
        // _OptimizedBEXY_ will fail due to an out of range error. In this situation,
        // releasing the accumulator qubit will throw an error as it will be in the One state.
    }

    // Subroutine of OptimizedBEXY.
    internal function _OptimizedBEXY_ (targetIndex : Int) : ((Qubit, Qubit, Qubit[]) => Unit is Adj + Ctl) {
        //Message($"OptimizedBEXY {targetIndex}");
        return _OptimizedBEXY(targetIndex, _, _, _);
    }


    // Subroutine of OptimizedBEXY.
    internal operation _OptimizedBEXY(targetIndex : Int, pauliBasis : Qubit, accumulator : Qubit, targetRegister : Qubit[]) : Unit is Adj {

        body (...) {
            // This should always be called as a controlled operation.
            fail "_OptimizedBEXY should always be called as a controlled operation.";
        }

        controlled (ctrl, ...) {
            if Length(targetRegister) <= targetIndex {
                fail "targetIndex out of range.";
            }

            Controlled X(ctrl, accumulator);
            ApplyWithCA(Controlled Adjoint S([pauliBasis], _), Controlled X(ctrl, _), targetRegister[targetIndex]);
            Controlled Z([accumulator], targetRegister[targetIndex]);
        }

    }


    /// # Summary
    /// Applies a Z operation to a qubit indicated by the state of another
    /// register.
    ///
    /// # Description
    /// The operation can be represented by a unitary matrix $U$ that applies
    /// the @"Microsoft.Quantum.Intrinsic.Z" operation on a qubit $p$
    /// conditioned on an index state $\ket{p}$. That is,
    /// $$
    /// \begin{align}
    ///     U\ket{p}\ket{\psi} = \ket{p}Z\_p\ket{\psi}.
    /// \end{align}
    /// $$
    ///
    /// # Input
    /// ## indexRegister
    /// A register in the state $\ket{p}$, determining the qubit on which $Z$ is applied.
    /// ## targetRegister
    /// Register of qubits on which the Pauli operators are applied.
    operation SelectZ (indexRegister : LittleEndian, targetRegister : Qubit[]) : Unit is Adj + Ctl {
        let unitaryGenerator = (Length(targetRegister), CurriedOpCA(ApplyToElementCA(Z, _, _)));
        MultiplexOperationsFromGenerator(unitaryGenerator, indexRegister, targetRegister);
        // If indexRegister encodes an integer that is larger than Length(targetRegister),
        // _SelectZ_ will fail due to an out of range error. In this situation,
        // releasing the accumulator qubit will throw an error as it will be in the One state.
    }

    internal operation _JordanWignerSelect_(
        signQubit : Qubit,
        selectZControlRegisters : Qubit[],
        OptimizedBEControlRegisters : Qubit[],
        pauliBases : Qubit[],
        indexRegisters : LittleEndian[],
        targetRegister : Qubit[]
    ) : Unit is Adj + Ctl {
        Z(signQubit);

        for idxRegister in IndexRange(OptimizedBEControlRegisters) {
            Controlled OptimizedBEXY([OptimizedBEControlRegisters[idxRegister]], (pauliBases[idxRegister], indexRegisters[idxRegister], targetRegister));
        }

        for idxRegister in IndexRange(selectZControlRegisters) {
            Controlled SelectZ([selectZControlRegisters[idxRegister]], (indexRegisters[idxRegister], targetRegister));
        }
    }


    internal function _JordanWignerSelectQubitCount_ (nZ : Int, nMaj : Int, nIdxRegQubits : Int) : (Int, (Int, Int, Int, Int, Int[])) {
        let signQubit = 1;
        let selectZControlRegisters = nZ;
        let OptimizedBEControlRegisters = nMaj;
        let pauliBases = nMaj;
        let indexRegisters = ConstantArray(Max([nMaj, nZ]), nIdxRegQubits);
        let nTotal = ((1 + nZ) + 2 * nMaj) + Max([nZ, nMaj]) * nIdxRegQubits;
        return (nTotal, (signQubit, selectZControlRegisters, OptimizedBEControlRegisters, pauliBases, indexRegisters));
    }


    internal function _JordanWignerSelectQubitManager_ (nZ : Int, nMaj : Int, nIdxRegQubits : Int, ctrlRegister : Qubit[], targetRegister : Qubit[]) : ((Qubit, Qubit[], Qubit[], Qubit[], LittleEndian[], Qubit[]), Qubit[]) {
        let (nTotal, (a, b, c, d, e)) = _JordanWignerSelectQubitCount_(nZ, nMaj, nIdxRegQubits);
        let split = [a, b, c, d] + e;
        let registers = Partitioned(split, ctrlRegister);
        let signQubit = registers[0];
        let selectZControlRegisters = registers[1];
        let OptimizedBEControlRegisters = registers[2];
        let pauliBases = registers[3];
        let indexRegistersTmp = registers[4 .. (4 + Length(e)) - 1];
        let rest = registers[Length(registers) - 1];
        mutable indexRegisters = new LittleEndian[Length(e)];

        for idx in IndexRange(e) {
            set indexRegisters w/= idx <- LittleEndian(indexRegistersTmp[idx]);
        }

        return ((signQubit[0], selectZControlRegisters, OptimizedBEControlRegisters, pauliBases, indexRegisters, targetRegister), rest);
    }

}
