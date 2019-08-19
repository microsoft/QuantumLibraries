// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Simulation {
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Preparation;
    open Microsoft.Quantum.Arithmetic;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Math;
    open Microsoft.Quantum.Convert;
    open Microsoft.Quantum.Arrays;

    /// # Summary
    /// Extracts the coefficient of a Pauli term described by a `GeneratorIndex`.
    ///
    /// # Input
    /// ## generatorIndex
    /// `GeneratorIndex` type that encodes a Pauli term.
    ///
    /// # Output
    /// The coefficient of the term described by a `GeneratorIndex`.
    function PauliCoefficientFromGenIdx(generatorIndex: GeneratorIndex) : Double {
        return Head(Snd(generatorIndex::IndexData));
    }

    /// # Summary
    /// Extracts the Pauli string and its qubit indices of a Pauli term described
    /// by a `GeneratorIndex`.
    ///
    /// # Input
    /// ## generatorIndex
    /// `GeneratorIndex` type that encodes a Pauli term.
    ///
    /// # Output
    /// The Pauli string of the term described by a `GeneratorIndex`, and
    /// indices to the qubits it acts on.
    function PauliStringFromGenIdx(generatorIndex: GeneratorIndex) : (Pauli[], Int[]) {
        let (idxPaulis, coeff) = generatorIndex::IndexData;
        return (IntArrayAsPauliArray(idxPaulis), generatorIndex::Subsystems);
    }

    /// # Summary
	/// Creates a block-encoding unitary for a Hamiltonian.
	/// 
    /// The Hamiltonian $H=\sum_{j}\alpha_j P_j$ is described by a 
    /// sum of Pauli terms $P_j$, each with real coefficient $\alpha_j$.
    ///
    /// # Input
    /// ## generatorSystem
    /// A `GeneratorSystem` that describes $H$ as a sum of Pauli terms
    ///
    /// # Output
    /// ## First parameter
    /// The one-norm of coefficients $\alpha=\sum_{j}|\alpha_j|$.
    /// ## Second parameter
    /// A `BlockEncodingReflection` unitary $U$ of the Hamiltonian $H$. As this unitary 
    /// satisfies $U^2 = I$, it is also a reflection.
    ///
    /// # Remarks
    /// This is obtained by preparing and unpreparing the state $\sum_{j}\sqrt{\alpha_j/\alpha}\ket{j}$,
    /// and constructing a multiply-controlled unitary 
    /// <xref:microsoft.quantum.preparation.statepreparationpositivecoefficients> and
    /// <xref:microsoft.quantum.canon.multiplexoperationsfromgenerator>.
    function PauliBlockEncoding(generatorSystem: GeneratorSystem) : (Double, BlockEncodingReflection) {
        let statePrepUnitary = StatePreparationPositiveCoefficients;
        let multiplexer = MultiplexerFromGenerator; 
        return PauliBlockEncodingImpl(generatorSystem, statePrepUnitary, multiplexer);
        
    }

    /// # Summary
	/// Creates a block-encoding unitary for a Hamiltonian.
	/// 
    /// The Hamiltonian $H=\sum_{j}\alpha_j P_j$ is described by a 
    /// sum of Pauli terms $P_j$, each with real coefficient $\alpha_j$.
    ///
    /// # Input
    /// ## generatorSystem
    /// A `GeneratorSystem` that describes $H$ as a sum of Pauli terms
    /// ## statePrepUnitary
    /// A unitary operation $P$ that prepares $P\ket{0}=\sum_{j}\sqrt{\alpha_j}\ket{j}$ given
    /// an array of coefficients $\{\sqrt{\alpha}_j\}$.
    /// ## statePrepUnitary
    /// A unitary operation $V$ that applies the unitary $V_j$ controlled on index $\ket{j}$,
    /// given a function $f: j\rightarrow V_j$.
    ///
    /// # Output
    /// ## First parameter
    /// The one-norm of coefficients $\alpha=\sum_{j}|\alpha_j|$.
    /// ## Second parameter
    /// A `BlockEncodingReflection` unitary $U$ of the Hamiltonian $U$. As this unitary 
    /// satisfies $U^2 = I$, it is also a reflection.
    ///
    /// # Remarks
    /// Example operations the prepare and unpreparing the state $\sum_{j}\sqrt{\alpha_j/\alpha}\ket{j}$,
    /// and construct a multiply-controlled unitary are
    /// <xref:microsoft.quantum.canon.statepreparationpositivecoefficients> and
    /// <xref:microsoft.quantum.canon.multiplexoperationsfromgenerator>.
    function PauliBlockEncodingImpl(
        generatorSystem: GeneratorSystem,
        statePrepUnitary: (Double[] -> (LittleEndian => Unit is Adj + Ctl)), 
        multiplexer: ((Int, (Int -> (Qubit[] => Unit is Adj + Ctl))) -> ((LittleEndian, Qubit[]) => Unit is Adj + Ctl))) : (Double, BlockEncodingReflection)
    {
        let op = _IdxToCoeff(_, generatorSystem::Term, PauliCoefficientFromGenIdx);
        let coefficients = Mapped(op, RangeAsIntArray(0..generatorSystem::NTerms - 1));
        let oneNorm = PowD(PNorm(2.0, coefficients), 2.0);
        let unitaryGenerator = (generatorSystem::NTerms, _IdxToUnitary(_, generatorSystem::Term, _PauliLCUUnitary));
        let statePreparation = statePrepUnitary(coefficients);
        let selector = multiplexer(unitaryGenerator);
        let blockEncoding = BlockEncodingReflection(BlockEncoding(_BlockEncodingFromBEandQubit(BlockEncodingByLCU(statePreparation, selector),_,_)));
        return (oneNorm, blockEncoding);
    }

    /// # Summary
    /// Used in implementation of `PauliBlockEncoding`
    /// # See Also
    /// - Microsoft.Quantum.Canon.PauliBlockEncoding
    function _IdxToCoeff(idx: Int, genFun: (Int -> GeneratorIndex), genIdxToCoeff: (GeneratorIndex -> Double)) : Double {
        return Sqrt(AbsD(genIdxToCoeff(genFun(idx))));
    }

    /// # Summary
    /// Used in implementation of `PauliBlockEncoding`
    /// # See Also
    /// - Microsoft.Quantum.Canon.PauliBlockEncoding
    function _IdxToUnitary(idx: Int, genFun: (Int -> GeneratorIndex), genIdxToUnitary: (GeneratorIndex -> (Qubit[] => Unit is Adj + Ctl))) : (Qubit[] => Unit is Adj + Ctl) {
        return genIdxToUnitary(genFun(idx));
    }


    /// # Summary
    /// Used in implementation of `PauliBlockEncoding`
    /// # See Also
    /// - Microsoft.Quantum.Canon.PauliBlockEncoding
    function _PauliLCUUnitary(generatorIndex: GeneratorIndex) : (Qubit[] => Unit is Adj + Ctl) {
        return _ApplyPauliLCUUnitary(generatorIndex,_);
    }

    /// # Summary
    /// Used in implementation of `PauliBlockEncoding`
    /// # See Also
    /// - Microsoft.Quantum.Canon.PauliBlockEncoding
    operation _ApplyPauliLCUUnitary(generatorIndex: GeneratorIndex, qubits: Qubit[]) : Unit
    is Adj + Ctl {
        let (idxPaulis, coeff) = generatorIndex::IndexData;
        let pauliString = IntArrayAsPauliArray(idxPaulis);
        let pauliQubits = Subarray(generatorIndex::Subsystems, qubits);

        ApplyPauli(pauliString, pauliQubits);

        if (coeff[0] < 0.0) {
            // -1 phase
            Exp([PauliI], PI(), [Head(pauliQubits)]);
        }
    }

}
