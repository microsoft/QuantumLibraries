// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Simulation {
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Arithmetic;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Math;

    // This library defines types for block-encoding an arbitrary operator.
    // General operations for performing block-encoding are provided,
    // together with operations for creating a quantum walk from This
    // block-encoding.

    /// # Summary
    /// Represents a unitary where an arbitrary operator of
    /// interest is encoded in the top-left block.
    ///
    /// # Description
    /// A unitary operation that can be represented by a unitary matrix $U$
    /// where an arbitrary operator $H$ of
    /// interest that acts on the system register `s` is encoded in the top-
    /// left block corresponding to auxiliary state $\ket{0}_a$. That is,
    ///
    /// $$
    /// \begin{align}
    /// (\bra{0}_a\otimes I_s)U(\ket{0}_a\otimes I_s) = H
    /// \end{align}
    /// $$.
    ///
    /// The inputs to this callable are:
    /// - An array of qubits representing the auxiliary register acted on by $U$.
    ///   Note that the action of $U$ is only defined when this register is
    ///   in the state $\ket{0}_a$.
    /// - An array of qubits representing the system register acted on by $H$.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Simulation.BlockEncodingReflection
    /// - Microsoft.Quantum.Simulation.TimeDependentBlockEncoding
    newtype BlockEncoding = ((Qubit[], Qubit[]) => Unit is Adj + Ctl);

    /// # Summary
    /// Represents a `BlockEncoding` that is also a reflection.
    ///
    /// # Description
    /// A unitary operation that can be represented by a unitary matrix $U$
    /// where an arbitrary operator $H$ of
    /// interest that acts on the system register `s` is encoded in the top-
    /// left block corresponding to auxiliary state $\ket{0}_a$. That is,
    ///
    /// $$
    /// \begin{align}
    /// (\bra{0}_a\otimes I_s)U(\ket{0}_a\otimes I_s) = H
    /// \end{align}
    /// $$.
    ///
    /// The inputs to this callable are:
    /// - An array of qubits representing the auxiliary register acted on by $U$.
    ///   Note that the action of $U$ is only defined when this register is
    ///   in the state $\ket{0}_a$.
    /// - An array of qubits representing the system register acted on by $H$.
    ///
    ///
    /// # See Also
    /// - Microsoft.Quantum.Simulation.BlockEncoding
    newtype BlockEncodingReflection = BlockEncoding;

    /// # Summary
    /// Represents a `BlockEncoding` that is controlled by a clock register.
    ///
    /// # Description
    /// An operation marked as a `TimeDependentBlockEncoding` can be
    /// represented by a unitary matrix $U$ controlled by a state
    /// $\ket{k}_d$ in clock register `d` such that an arbitrary operator $H_k$ of
    /// interest that acts on the system register `s` is encoded in the top-
    /// left block corresponding to auxiliary state $\ket{0}_a$. That is,
    ///
    /// $$
    /// \begin{align}
    /// (\bra{0}\_a\otimes I_{ds})U(\ket{0}\_a\otimes I_{ds}) = \sum_{k}\ket{k}\bra{k}\_d\otimes H_k.
    /// \end{align}
    /// $$.
    ///
    /// The inputs to the operation wrapped by this user-defined type are:
    /// - An array of qubits representing the time register that controls $H_k$.
    /// - An array of qubits representing the auxiliary register acted on by $U$.
    ///   The action of $U$ is only defined when this is $\ket{0}_a$.
    /// - An array of qubits representing the system register acted on by $H$.
    newtype TimeDependentBlockEncoding = ((Qubit[], Qubit[], Qubit[]) => Unit is Adj + Ctl);

    /// # Summary
    /// Converts a `BlockEncoding` into an equivalent `BLockEncodingReflection`.
    ///
    /// That is, given a `BlockEncoding` unitary $U$ that encodes some
    /// operator $H$ of interest, converts it into a `BlockEncodingReflection` $U'$ that
    /// encodes the same operator, but also satisfies $U'^\dagger = U'$.
    /// This increases the size of the auxiliary register of $U$ by one qubit.
    ///
    /// # Input
    /// ## blockEncoding
    /// An operation to be converted into a reflection, and that is represented
    /// by a unitary matrix $U$.
    ///
    /// # Output
    /// A operation represented by a unitary matrix $U'$ acting jointly on
    /// registers `a` and `s` that block-encodes $H$, and
    /// satisfies $U'^\dagger = U'$.
    ///
    /// # Remarks
    /// This increases the size of the auxiliary register of $U$ by one qubit.
    ///
    /// # References
    /// - Hamiltonian Simulation by Qubitization
    ///   Guang Hao Low, Isaac L. Chuang
    ///   https://arxiv.org/abs/1610.06546
    ///
    /// # See Also
    /// - Microsoft.Quantum.Simulation.BlockEncoding
    /// - Microsoft.Quantum.Simulation.BlockEncodingReflection
    function BlockEncodingToReflection(blockEncoding : BlockEncoding)
    : BlockEncodingReflection {
        return BlockEncodingReflection(BlockEncoding(ApplyBlockEncodingAsReflection(blockEncoding, _, _)));
    }

    /// # Summary
    /// Implementation of `BlockEncodingToReflection`.
    internal operation ApplyBlockEncodingAsReflection(blockEncoding: BlockEncoding, auxiliary: Qubit[], system: Qubit[])
    : Unit is Adj + Ctl {
        let prep = auxiliary[0];
        let blockEncodingAux = Rest(auxiliary);
        let op1 = Controlled blockEncoding!(_, (blockEncodingAux, system));
        let op0 = ApplyToEachCA(H,_);
        ApplyWithCA(op0, ApplyWithCA(op1, ApplyToEachCA(X,_), _), [prep]);
    }

    /// # Summary
    /// Converts a block-encoding reflection into a quantum walk.
    ///
    /// # Description
    /// Given a `BlockEncodingReflection` represented by a unitary $U$
    /// that encodes some operator $H$ of interest, converts it into a quantum
    /// walk $W$ containing the spectrum of $\pm e^{\pm i\sin^{-1}(H)}$.
    ///
    /// # Input
    /// ## blockEncoding
    /// A `BlockEncodingReflection` unitary $U$ to be converted into a Quantum
    /// walk.
    ///
    /// # Output
    /// A quantum walk $W$ acting jointly on registers `a` and `s` that block-
    /// encodes $H$, and contains the spectrum of $\pm e^{\pm i\sin^{-1}(H)}$.
    ///
    /// # References
    /// - Hamiltonian Simulation by Qubitization
    ///   Guang Hao Low, Isaac L. Chuang
    ///   https://arxiv.org/abs/1610.06546
    ///
    /// # See Also
    /// - Microsoft.Quantum.Simulation.BlockEncoding
    /// - Microsoft.Quantum.Simulation.BlockEncodingReflection
    function QuantumWalkByQubitization(blockEncoding : BlockEncodingReflection)
    : ((Qubit[], Qubit[]) => Unit is Adj + Ctl) {
        return ApplyQuantumWalkByQubitization(blockEncoding, _, _);
    }

    /// # Summary
    /// Implementation of `Qubitization`.
    internal operation ApplyQuantumWalkByQubitization(
        blockEncoding : BlockEncodingReflection,
        auxiliary : Qubit[],
        system : Qubit[]
    )
    : Unit is Adj + Ctl {
        Exp([PauliI], -0.5 * PI(), [Head(system)]);
        RAll0(PI(), auxiliary);
        // NB: We unwrap twice here, since BlockEncodingReflection wraps BlockEncoding.
        blockEncoding!!(auxiliary, system);
    }

    /// # Summary
    /// Encodes an operator of interest into a `BlockEncoding`.
    ///
    /// This constructs a `BlockEncoding` unitary $U=P\cdot V\cdot P^\dagger$ that encodes some
    /// operator $H = \sum_{j}|\alpha_j|U_j$ of interest that is a linear combination of
    /// unitaries. Typically, $P$ is a state preparation unitary such that
    /// $P\ket{0}\_a=\sum_j\sqrt{\alpha_j/\|\vec\alpha\|\_2}\ket{j}\_a$,
    /// and $V=\sum_{j}\ket{j}\bra{j}\_a\otimes U_j$.
    ///
    /// # Input
    /// ## statePreparation
    /// A unitary $P$ that prepares some target state.
    /// ## selector
    /// A unitary $V$ that encodes the component unitaries of $H$.
    ///
    /// # Output
    /// A unitary $U$ acting jointly on registers `a` and `s` that block-
    /// encodes $H$, and satisfies $U^\dagger = U$.
    ///
    /// # Remarks
    /// This `BlockEncoding` implementation gives it the properties of a
    /// `BlockEncodingReflection`.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Simulation.BlockEncoding
    /// - Microsoft.Quantum.Simulation.BlockEncodingReflection
    function BlockEncodingByLCU<'T,'S>(
        statePreparation: ('T => Unit is Adj + Ctl),
        selector: (('T, 'S) => Unit is Adj + Ctl)
    )
    : (('T, 'S) => Unit is Adj + Ctl) {
        return ApplyBlockEncodingByLCU(statePreparation, selector, _, _);
    }

    /// # Summary
    /// Implementation of `BlockEncodingByLCU`.
    internal operation ApplyBlockEncodingByLCU<'T,'S>(
        statePreparation: ('T => Unit is Adj + Ctl),
        selector: (('T, 'S) => Unit is Adj + Ctl),
        auxiliary: 'T,
        system: 'S
    )
    : Unit is Adj + Ctl {
        ApplyWithCA(statePreparation, selector(_, system), auxiliary);
    }

    /// # Summary
    /// Encodes an operator of interest into a `BlockEncodingReflection`.
    ///
    /// This constructs a `BlockEncodingReflection` unitary $U=P\cdot V\cdot P^\dagger$ that encodes some
    /// operator $H=\sum_{j}|\alpha_j|U_j$ of interest that is a linear combination of
    /// unitaries. Typically, $P$ is a state preparation unitary such that
    /// $P\ket{0}\_a\sum_j\sqrt{\alpha_j/\|\vec\alpha\|\_2}\ket{j}\_a$,
    /// and $V=\sum_{j}\ket{j}\bra{j}\_a\otimes U_j$.
    ///
    /// # Input
    /// ## statePreparation
    /// A unitary $P$ that prepares some target state.
    /// ## selector
    /// A unitary $V$ that encodes the component unitaries of $H$.
    ///
    /// # Output
    /// A unitary $U$ acting jointly on registers `a` and `s` that block-
    /// encodes $H$, and satisfies $U'^{-1} = U$.
    ///
    /// # Remarks
    /// This `BlockEncoding` implementation gives it the properties of a
    /// `BlockEncodingReflection`.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Simulation.BlockEncoding
    /// - Microsoft.Quantum.Simulation.BlockEncodingReflection
    function BlockEncodingReflectionByLCU(
        statePreparation : (Qubit[] => Unit is Adj + Ctl),
        selector : ((Qubit[], Qubit[]) => Unit is Adj + Ctl)
    )
    : BlockEncodingReflection {
        return BlockEncodingToReflection(BlockEncoding(BlockEncodingByLCU(statePreparation, selector)));
    }


    /// # Summary
    /// Conversion of ((LittleEndian, Qubit[]) => () is Adj + Ctl) to BlockEncoding
    internal operation ApplyBlockEncodingFromBEandQubit(
        op: ((LittleEndian, Qubit[]) => Unit is Adj + Ctl),
        auxiliary: Qubit[],
        system: Qubit[]
    )
    : Unit is Adj + Ctl {
        op(LittleEndian(auxiliary), system);
    }

}
