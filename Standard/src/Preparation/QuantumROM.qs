// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Preparation {
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Arithmetic;
    open Microsoft.Quantum.Convert;
    open Microsoft.Quantum.Math;
    open Microsoft.Quantum.Arrays;

    /// # Summary
    /// Returns an operation that prepares a given mixed state.
    ///
    /// # Description
    /// Uses the Quantum ROM technique to represent a given density matrix,
    /// returning that representation as a state preparation operation.
    ///
    /// In particular, given a list of $N$ coefficients $\alpha_j$, this
    /// function returns an operation that uses the Quantum ROM technique to
    /// prepare an approximation
    /// $$
    /// \begin{align}
    ///     \tilde\rho = \sum_{j = 0}^{N - 1} p_j \ket{j}\bra{j}
    /// \end{align}
    /// $$
    /// of the mixed state
    /// $$
    /// \begin{align}
    ///     \rho = \sum_{j = 0}^{N-1}\ frac{|alpha_j|}{\sum_k |\alpha_k|} \ket{j}\bra{j},
    /// \end{align}
    /// $$
    /// where each $p_j$ is an approximation to the given coefficient $\alpha_j$
    /// such that
    /// $$
    /// \begin{align}
    ///     \left| p_j - \frac{ |\alpha_j| }{ \sum_k |\alpha_k| } \le \frac{\epsilon}{N}
    /// \end{align}
    /// $$
    /// for each $j$.
    ///
    /// When passed an index register and a register of garbage qubits,
    /// initially in the state $\ket{0} \ket{00\cdots 0}, the returned operation
    /// prepares both registers into the purification of $\tilde \rho$,
    /// $$
    /// \begin{align}
    ///     \sum_{j=0}^{N-1} \sqrt{p_j} \ket{j}\ket{\text{garbage}_j},
    /// \end{align}
    /// $$
    /// such that resetting and deallocating the garbage register enacts the
    /// desired preparation to within the target error $\epsilon$.
    ///
    /// # Input
    /// ## targetError
    /// The target error $\epsilon$.
    /// ## coefficients
    /// Array of $N$ coefficients specifying the probability of basis states.
    /// Negative numbers $-\alpha_j$ will be treated as positive $|\alpha_j|$.
    ///
    /// # Output
    /// An operation that prepares $\tilde \rho$ as a purification onto a joint
    /// index and garbage register.
    ///
    /// # Remarks
    /// The coefficients provided to this operation are normalized following the
    /// 1-norm, such that the coefficients are always considered to describe a
    /// valid categorical probability distribution.
    ///
    /// # Example
    /// The following code snippet prepares an purification of the $3$-qubit state
    /// $\rho=\sum_{j=0}^{4}\frac{|alpha_j|}{\sum_k |\alpha_k|}\ket{j}\bra{j}$, where
    /// $\vec\alpha=(1.0, 2.0, 3.0, 4.0, 5.0)$, and the target error is
    /// $10^{-3}$:
    /// ```Q#
    /// let coefficients = [1.0, 2.0, 3.0, 4.0, 5.0];
    /// let targetError = 1e-3;
    /// let preparation = MixedStatePreparation(targetError, coefficients);
    /// using (indexRegister = Qubit[preparation::Requirements::NIndexQubits]) {
    ///     using (garbageRegister = Qubit[preparation::Requirements::NGarbageQubits]) {
    ///         preparation::Prepare(LittleEndian(indexRegister), garbageRegister);
    ///     }
    /// }
    /// ```
    ///
    /// # References
    /// - Encoding Electronic Spectra in Quantum Circuits with Linear T Complexity
    ///   Ryan Babbush, Craig Gidney, Dominic W. Berry, Nathan Wiebe, Jarrod McClean, Alexandru Paler, Austin Fowler, Hartmut Neven
    ///   https://arxiv.org/abs/1805.03662
    function MixedStatePreparation(targetError : Double, coefficients : Double[])
    : MixedPreparationOperation {
        let nBitsPrecision = -Ceiling(Lg(0.5 * targetError)) + 1;
        let (oneNorm, keepCoeff, altIndex) = _QuantumROMDiscretization(nBitsPrecision, coefficients);
        let nCoeffs = Length(coefficients);
        let nBitsIndices = Ceiling(Lg(IntAsDouble(nCoeffs)));

        let op =  _QuantumROMImpl(nBitsPrecision, nCoeffs, nBitsIndices, keepCoeff, altIndex, _, _);
        let qubitCounts = MixedStatePreparationRequirements(targetError, nCoeffs);
        return MixedPreparationOperation(qubitCounts, oneNorm, op);
    }

    /// # Summary
    /// Returns the total number of qubits that must be allocated
    /// to the operation returned by
    /// @"microsoft.quantum.preparation.mixedstatepreparation".
    ///
    /// # Input
    /// ## targetError
    /// The target error $\epsilon$.
    /// ## nCoefficients
    /// The number of coefficients to be specified in preparing a mixed state.
    ///
    /// # Output
    /// A description of how many qubits are required in total, and for each of
    /// the index and garbage registers used by the
    /// @"microsoft.quantum.preparation.mixedstatepreparation" function.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Preparation.MixedStatePreparation
    function MixedStatePreparationRequirements(targetError : Double, nCoefficients : Int)
    : MixedPreparationRequirements {
        let nBitsPrecision = -Ceiling(Lg(0.5*targetError)) + 1;
        let nIndexQubits = Ceiling(Lg(IntAsDouble(nCoefficients)));
        let nGarbageQubits = nIndexQubits + 2 * nBitsPrecision + 1;
        let nTotal = nGarbageQubits + nIndexQubits;
        return MixedPreparationRequirements(nTotal, (nIndexQubits, nGarbageQubits));
    }

    // Implementation step of `QuantumROM`. This splits a single
    // qubit array into the subarrays required by the operation.
    function _PartitionedForQuantumROM(targetError: Double, nCoeffs: Int, qubits: Qubit[])
    : ((LittleEndian, Qubit[]), Qubit[]) {
        let requirements = MixedStatePreparationRequirements(targetError, nCoeffs);
        let registers = Partitioned(
            [requirements::NIndexQubits, requirements::NGarbageQubits],
            qubits
        );
        return ((LittleEndian(registers[0]), registers[1]), registers[2]);
    }

    // Classical processing
    // This discretizes the coefficients such that
    // |coefficient[i] * oneNorm - discretizedCoefficient[i] * discreizedOneNorm| * nCoeffs <= 2^{1-bitsPrecision}.
    function _QuantumROMDiscretization(bitsPrecision: Int, coefficients: Double[])
    : (Double, Int[], Int[]) {
        let oneNorm = PNorm(1.0, coefficients);
        let nCoefficients = Length(coefficients);
        if (bitsPrecision > 31) {
            fail $"Bits of precision {bitsPrecision} unsupported. Max is 31.";
        }
        if (nCoefficients <= 1) {
            fail $"Cannot prepare state with less than 2 coefficients.";
        }
        if (oneNorm == 0.0) {
            fail $"State must have at least one coefficient > 0";
        }

        let barHeight = 2 ^ bitsPrecision - 1;

        mutable altIndex = RangeAsIntArray(0..nCoefficients - 1);
        mutable keepCoeff = Mapped(
            _QuantumROMDiscretizationRoundCoefficients(_, oneNorm, nCoefficients, barHeight),
            coefficients
        );

        // Calculate difference between number of discretized bars vs. maximum
        mutable bars = 0;
        for (idxCoeff in IndexRange(keepCoeff)) {
            set bars += keepCoeff[idxCoeff] - barHeight;
        }
        //Message($"Excess bars {bars}.");
        // Uniformly distribute excess bars across coefficients.
        for (idx in 0..AbsI(bars) - 1) {
            set keepCoeff w/= idx <- keepCoeff[idx] + bars > 0 ? 1 | -1;
        }

        mutable barSink = new Int[0];
        mutable barSource = new Int[0];

        for (idxCoeff in IndexRange(keepCoeff)) {
            if (keepCoeff[idxCoeff] > barHeight) {
                set barSource += [idxCoeff];
            } elif (keepCoeff[idxCoeff] < barHeight) {
                set barSink += [idxCoeff];
            }
        }

        for (rep in 0..nCoefficients * 10) {
            if (Length(barSink) > 0 and Length(barSource) > 0) {
                let idxSink = Tail(barSink);
                let idxSource = Tail(barSource);
                set barSink = Most(barSink);
                set barSource = Most(barSource);

                set keepCoeff w/= idxSource <- keepCoeff[idxSource] - barHeight + keepCoeff[idxSink];
                set altIndex w/= idxSink <- idxSource;

                if (keepCoeff[idxSource] < barHeight) {
                    set barSink += [idxSource];
                } elif (keepCoeff[idxSource] > barHeight) {
                    set barSource += [idxSource];
                }
            } elif (Length(barSource) > 0) {
                //Message($"rep: {rep}, nBarSource {nBarSource}.");
                let idxSource = Tail(barSource);
                set barSource = Most(barSource);
                set keepCoeff w/= idxSource <- barHeight;
            } else {
                return (oneNorm, keepCoeff, altIndex);
            }
        }

        return (oneNorm, keepCoeff, altIndex);
    }

    // Used in QuantumROM implementation.
    function _QuantumROMDiscretizationRoundCoefficients(coefficient: Double, oneNorm: Double, nCoefficients: Int, barHeight: Int) : Int {
        return Round((AbsD(coefficient) / oneNorm) * IntAsDouble(nCoefficients) * IntAsDouble(barHeight));
    }

    // Used in QuantumROM implementation.
    operation _QuantumROMImpl(
        nBitsPrecision: Int, nCoeffs: Int, nBitsIndices: Int, keepCoeff: Int[],
        altIndex: Int[], indexRegister: LittleEndian, garbageRegister: Qubit[]
    )
    : Unit is Adj + Ctl {
        let unitaryGenerator = (nCoeffs, _QuantumROMWriteBitStringUnitary(_, keepCoeff, altIndex));
        let partitionedRegisters = Partitioned([nBitsIndices, nBitsPrecision, nBitsPrecision, 1], garbageRegister);

        let altIndexRegister = LittleEndian(partitionedRegisters[0]);
        let keepCoeffRegister = LittleEndian(partitionedRegisters[1]);
        let uniformKeepCoeffRegister = LittleEndian(partitionedRegisters[2]);
        let flagQubit = partitionedRegisters[3][0];

        // Create uniform superposition over index and alt coeff register.
        PrepareUniformSuperposition(nCoeffs, indexRegister);
        ApplyToEachCA(H, uniformKeepCoeffRegister!);

        // Write bitstrings to altIndex and keepCoeff register.
        MultiplexOperationsFromGenerator(unitaryGenerator, indexRegister, (keepCoeffRegister, altIndexRegister));

        // Perform comparison
        CompareUsingRippleCarry(uniformKeepCoeffRegister, keepCoeffRegister, flagQubit);

        // Swap in register based on comparison
        Controlled ApplyToEachCA(
            [flagQubit],
            (
                SWAP,
                Zip(Reversed(indexRegister!), altIndexRegister!)
            )
        );
        // for (idx in 0..nBitsIndices - 1) {
        //     (Controlled SWAP)([flagQubit], (indexRegister![nBitsIndices - idx - 1], altIndexRegister![idx]));
        // }
    }
    // Used in QuantumROM implementation.
    function _QuantumROMWriteBitStringUnitary(idx : Int, keepCoeff : Int[], altIndex : Int[])
    : ((LittleEndian, LittleEndian) => Unit is Adj + Ctl) {
        return _QuantumROMWriteBitString(idx, keepCoeff, altIndex, _, _);
    }

    // Used in QuantumROM implementation.
    operation _QuantumROMWriteBitString(
        idx: Int, keepCoeff: Int[], altIndex: Int[],
        keepCoeffRegister: LittleEndian, altIndexRegister: LittleEndian
    )
    : Unit is Adj + Ctl {
        ApplyXorInPlace(keepCoeff[idx], keepCoeffRegister);
        ApplyXorInPlace(altIndex[idx], altIndexRegister);
    }



}
