// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Simulation {
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Arithmetic;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Arrays;

    // This library defines types for block-encoding an arbitrary operator.
    // General operations for performing block-encoding are provided,
    // together with operations for creating a quantum walk from This
    // block-encoding.

    /// # Summary
    /// Represents a unitary where an arbitrary operator of 
    /// interest is encoded in the top-left block.
	///
    /// That is, a `BlockEncoding` is a unitary $U$ where an arbitrary operator $H$ of 
    /// interest that acts on the system register `s` is encoded in the top-
    /// left block corresponding to auxiliary state $\ket{0}_a$. That is,
    ///
    /// $$
    /// \begin{align}
    /// (\bra{0}_a\otimes I_s)U(\ket{0}_a\otimes I_s) = H
    /// \end{align}
    /// $$.
    ///
    /// # Input
    /// ## First Parameter
    /// An array of qubits representing the auxiliary register acted on by $U$.
    /// The action of $U$ is only defined when this is $\ket{0}_a$.
    /// ## Second Parameter
    /// An array of qubits representing the system register acted on by $H$.
    ///
    /// # Output
    /// A unitary $U$ acting jointly on registers `a` and `s` that block-
    /// encodes $H$.
    newtype BlockEncoding = ((Qubit[], Qubit[]) => Unit : Adjoint, Controlled);

    /// # Summary
    /// Represents a `BlockEncoding` that is also a reflection.
    ///
    /// # Input
    /// ## First Parameter
    /// An array of qubits representing the auxiliary register acted on by $U$.
    /// The action of $U$ is only defined when this is $\ket{0}_a$.
    /// ## Second Parameter
    /// An array of qubits representing the system register acted on by $H$.
    ///
    /// # Output
    /// A unitary $U$ acting jointly on registers `a` and `s` that block-
    /// encodes $H$.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.BlockEncoding
    newtype BlockEncodingReflection = BlockEncoding;

    /// # Summary
	/// Represents a `BlockEncoding` that is controlled by a clock register.
	/// 
    /// That is, a `TimeDependentBlockEncoding` is a unitary $U$ controlled by a state 
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
    /// # Input
    /// ## First Parameter
    /// An array of qubits representing the time register that controls $H_k$.
    /// ## Second Parameter
    /// An array of qubits representing the auxiliary register acted on by $U$.
    /// The action of $U$ is only defined when this is $\ket{0}_a$.
    /// ## Third Parameter
    /// An array of qubits representing the system register acted on by $H$.
    ///
    /// # Output
    /// A unitary $U$ acting jointly on registers `d`, `a`, and `s` that block-
    /// encodes $H_k$.
    newtype TimeDependentBlockEncoding = ((Qubit[], Qubit[], Qubit[]) => Unit : Adjoint, Controlled);

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
    /// A `BlockEncoding` unitary $U$ to be converted into a reflection.
    ///
    /// # Output
    /// A unitary $U'$ acting jointly on registers `a` and `s` that block-
    /// encodes $H$, and satisfies $U'^\dagger = U'$.
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
    /// - Microsoft.Quantum.Canon.BlockEncoding
    /// - Microsoft.Quantum.Canon.BlockEncodingReflection
    function BlockEncodingToReflection(blockEncoding: BlockEncoding) : BlockEncodingReflection
    {
        return BlockEncodingReflection(BlockEncoding(BlockEncodingToReflection_(blockEncoding, _, _)));
    }
    
    /// # Summary
    /// Implementation of `BlockEncodingToReflection`.
    operation BlockEncodingToReflection_(blockEncoding: BlockEncoding, auxiliary: Qubit[], system: Qubit[]) : Unit {
        body (...) {
            let prep = auxiliary[0];
            let blockEncodingAncilla = Rest(auxiliary);
            let op1 = Controlled blockEncoding!(_, (blockEncodingAncilla, system));
            let op0 = ApplyToEachCA(H,_);
            ApplyWithCA(op0, ApplyWithCA(op1, ApplyToEachCA(X,_), _), [prep]);
        }
        adjoint auto;
        controlled auto;
        adjoint controlled auto;
    }

    /// # Summary
	/// Converts a `BlockEncodingReflection` into a quantum walk.
	/// 
    /// That is, given a `BlockEncodingReflection` unitary $U$ 
    /// that encodes some operator $H$ of interest, converts it into a quantum walk
    /// $W$ containing the spectrum of $\pm e^{\pm i\sin^{-1}(H)}$.
    ///
    /// # Input
    /// ## blockEncodingReflection
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
    /// - Microsoft.Quantum.Canon.BlockEncoding
    /// - Microsoft.Quantum.Canon.BlockEncodingReflection
    function QuantumWalkByQubitization(blockEncoding: BlockEncodingReflection) : ((Qubit[], Qubit[]) => Unit : Adjoint, Controlled) {
        return QuantumWalkByQubitization_(blockEncoding, _, _);
    }

    /// # Summary
    /// Implementation of `Qubitization`.
    operation QuantumWalkByQubitization_(blockEncoding: BlockEncodingReflection, auxiliary: Qubit[], system: Qubit[]) : Unit {    
        body (...){
            Exp([PauliI], -0.5 * Microsoft.Quantum.Extensions.Math.PI(), [system[0]]);
            RAll0(Microsoft.Quantum.Extensions.Math.PI(), auxiliary);
            // NB: We unwrap twice here, since BlockEncodingReflection wraps BlockEncoding.
            blockEncoding!!(auxiliary, system);    
        }
        adjoint auto;
        controlled auto;
        adjoint controlled auto;
    }

    /// # Summary
	/// Encodes an operator of interest into a `BlockEncoding`.
	/// 
    /// This constructs a `BlockEncoding` unitary $U=P\cdot V\cdot P^\dagger$ that encodes some
    /// operator $H=\sum_{j}|\alpha_j|U_j$ of interest that is a linear combination of
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
    /// - Microsoft.Quantum.Canon.BlockEncoding
    /// - Microsoft.Quantum.Canon.BlockEncodingReflection
    function BlockEncodingByLCU<'T,'S>(
        statePreparation: ('T => Unit : Adjoint, Controlled), 
        selector: (('T, 'S) => Unit : Adjoint, Controlled))
        : (('T, 'S) => Unit : Adjoint, Controlled) {
        return BlockEncodingByLCU_(statePreparation, selector, _, _);
    } 

    /// # Summary
    /// Implementation of `BlockEncodingByLCU`.
    operation BlockEncodingByLCU_<'T,'S>(
        statePreparation: ('T => Unit : Adjoint, Controlled), 
        selector: (('T, 'S) => Unit : Adjoint, Controlled),
        auxiliary: 'T, 
        system: 'S) 
        : Unit {
        body (...){
            ApplyWithCA(statePreparation, selector(_, system), auxiliary);
        }
        adjoint auto;
        controlled auto;
        adjoint controlled auto;
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
    /// - Microsoft.Quantum.Canon.BlockEncoding
    /// - Microsoft.Quantum.Canon.BlockEncodingReflection
    function BlockEncodingReflectionByLCU(
        statePreparation: (Qubit[] => Unit : Adjoint, Controlled), 
        selector: ((Qubit[], Qubit[]) => Unit : Adjoint, Controlled)
        ) : BlockEncodingReflection {
        return BlockEncodingToReflection(BlockEncoding(BlockEncodingByLCU(statePreparation, selector)));
    }


    /// # Summary
    /// Conversion of ((BigEndian, Qubit[]) => () : Adjoint, Controlled) to BlockEncoding
    operation BlockEncodingFromBEandQubit_(
        op: ((BigEndian, Qubit[]) => Unit : Adjoint, Controlled), 
        auxiliary: Qubit[], 
        system: Qubit[])  : Unit
    {
        body (...) {
            op(BigEndian(auxiliary), system);
        }
        adjoint auto;
        controlled auto;
        adjoint controlled auto;
    }
        
}
