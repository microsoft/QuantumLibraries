// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Preparation {
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Arithmetic;
    open Microsoft.Quantum.Convert;
    open Microsoft.Quantum.Math;
    open Microsoft.Quantum.Arrays;

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
    function QuantumROM(targetError: Double, coefficients: Double[])
    : (MixedStatePreparationRequirements, Double, ((LittleEndian, Qubit[]) => Unit is Adj + Ctl)) {
        let nBitsPrecision = -Ceiling(Lg(0.5 * targetError)) + 1;
        let positiveCoefficients = Mapped(AbsD, coefficients);
        let (oneNorm, keepCoeff, altIndex) = _QuantumROMDiscretization(nBitsPrecision, positiveCoefficients);
        let nCoeffs = Length(positiveCoefficients);
        let nBitsIndices = Ceiling(Lg(IntAsDouble(nCoeffs)));

        let op = PrepareQuantumROMStateWithoutSign(nBitsPrecision, nCoeffs, nBitsIndices, keepCoeff, altIndex, _, _);
        let qubitCounts = PurifiedMixedStateRequirements(targetError, nCoeffs);
        return (qubitCounts, oneNorm, op);
    }

    internal function SplitSign(coefficient : Double) : (Double, Int) {
        return (AbsD(coefficient), coefficient < 0.0 ? 1 | 0);
    }

    function QuantumROMWithSign(targetError : Double, coefficients : Double[])
    : (MixedStatePreparationRequirements, Double, ((LittleEndian, Qubit, Qubit[]) => Unit is Adj + Ctl)) {
        let nBitsPrecision = -Ceiling(Lg(0.5 * targetError)) + 1;
        let (positiveCoefficients, signs) = Unzipped(Mapped(SplitSign, coefficients));
        let (oneNorm, keepCoeff, altIndex) = _QuantumROMDiscretization(nBitsPrecision, positiveCoefficients);
        let nCoeffs = Length(positiveCoefficients);
        let nBitsIndices = Ceiling(Lg(IntAsDouble(nCoeffs)));

        let op = PrepareQuantumROMStateWithSign(nBitsPrecision, nCoeffs, nBitsIndices, keepCoeff, altIndex, signs, _, _, _);
        let qubitCounts = PurifiedMixedStateRequirements(targetError, nCoeffs);
        return (qubitCounts w/ NGarbageQubits <- qubitCounts::NGarbageQubits + 1, oneNorm, op);
    }

    /// # Summary
    /// Returns the total number of qubits that must be allocated
    /// in order to apply the operation returned by
    /// @"microsoft.quantum.preparation.purifiedmixedstate".
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
    /// @"microsoft.quantum.preparation.purifiedmixedstate" function.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Preparation.PurifiedMixedState
    function PurifiedMixedStateRequirements(targetError : Double, nCoefficients : Int)
    : MixedStatePreparationRequirements {
        let nBitsPrecision = -Ceiling(Lg(0.5*targetError)) + 1;
        let nIndexQubits = Ceiling(Lg(IntAsDouble(nCoefficients)));
        let nGarbageQubits = nIndexQubits + 2 * nBitsPrecision + 1;
        let nTotal = nGarbageQubits + nIndexQubits;
        return MixedStatePreparationRequirements(nTotal, (nIndexQubits, nGarbageQubits));
    }

    // Classical processing
    // This discretizes the coefficients such that
    // |coefficient[i] * oneNorm - discretizedCoefficient[i] * discretizedOneNorm| * nCoeffs <= 2^{1-bitsPrecision}.
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
            QuantumROMDiscretizationRoundCoefficients(_, oneNorm, nCoefficients, barHeight),
            coefficients
        );

        // Calculate difference between number of discretized bars vs. maximum
        mutable bars = 0;
        for (idxCoeff in IndexRange(keepCoeff)) {
            set bars += keepCoeff[idxCoeff] - barHeight;
        }

        // Uniformly distribute excess bars across coefficients.
        for (idx in 0..AbsI(bars) - 1) {
            set keepCoeff w/= idx <- keepCoeff[idx] + (bars > 0 ? -1 | +1);
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
    internal function QuantumROMDiscretizationRoundCoefficients(coefficient: Double, oneNorm: Double, nCoefficients: Int, barHeight: Int) : Int {
        return Round((AbsD(coefficient) / oneNorm) * IntAsDouble(nCoefficients) * IntAsDouble(barHeight));
    }

    // Used in QuantumROM implementation.
    internal function RoundedDiscretizationCoefficients(coefficient: Double, oneNorm: Double, nCoefficients: Int, barHeight: Int)
    : Int {
        return Round((AbsD(coefficient) / oneNorm) * IntAsDouble(nCoefficients) * IntAsDouble(barHeight));
    }

    // Used in QuantumROM implementation.
    internal operation PrepareQuantumROMState(nBitsPrecision: Int, nCoeffs: Int, nBitsIndices: Int, keepCoeff: Int[], altIndex: Int[], signs : Int[], indexRegister: LittleEndian, dataQubits : Qubit[], garbageRegister: Qubit[])
    : Unit is Adj + Ctl {
        let garbageIdx0 = nBitsIndices;
        let garbageIdx1 = garbageIdx0 + nBitsPrecision;
        let garbageIdx2 = garbageIdx1 + nBitsPrecision;
        let garbageIdx3 = garbageIdx2 + 1;

        let altIndexRegister = LittleEndian(garbageRegister[0..garbageIdx0 - 1]);
        let keepCoeffRegister = LittleEndian(garbageRegister[garbageIdx0..garbageIdx1 - 1]);
        let uniformKeepCoeffRegister = LittleEndian(garbageRegister[garbageIdx1..garbageIdx2 - 1]);
        let flagQubit = garbageRegister[garbageIdx3 - 1];
        let dataRegister = LittleEndian(dataQubits);
        let altDataRegister = LittleEndian(garbageRegister[garbageIdx3...]);

        // Create uniform superposition over index and alt coeff register.
        PrepareUniformSuperposition(nCoeffs, indexRegister);
        ApplyToEachCA(H, uniformKeepCoeffRegister!);

        // Write bitstrings to altIndex and keepCoeff register.
        let unitaryGenerator = (nCoeffs, QuantumROMBitStringWriterByIndex(_, keepCoeff, altIndex, signs));
        MultiplexOperationsFromGenerator(unitaryGenerator, indexRegister, (keepCoeffRegister, altIndexRegister, dataRegister, altDataRegister));

        // Perform comparison
        CompareUsingRippleCarry(uniformKeepCoeffRegister, keepCoeffRegister, flagQubit);

        let indexRegisterSize = Length(indexRegister!);

        // Swap in register based on comparison
        ApplyToEachCA((Controlled SWAP)([flagQubit], _), Zip(indexRegister! + dataRegister!, altIndexRegister! + altDataRegister!));
    }

    // # Remark
    // Application case for Maybe UDT
    internal operation PrepareQuantumROMStateWithoutSign(nBitsPrecision: Int, nCoeffs: Int, nBitsIndices: Int, keepCoeff: Int[], altIndex: Int[], indexRegister: LittleEndian, garbageRegister: Qubit[])
    : Unit is Adj + Ctl {
        PrepareQuantumROMState(nBitsPrecision, nCoeffs, nBitsIndices, keepCoeff, altIndex, new Int[0], indexRegister, new Qubit[0], garbageRegister);
    }

    // # Remark
    // Application case for Maybe UDT
    internal operation PrepareQuantumROMStateWithSign(nBitsPrecision: Int, nCoeffs: Int, nBitsIndices: Int, keepCoeff: Int[], altIndex: Int[], signs : Int[], indexRegister: LittleEndian, signQubit : Qubit, garbageRegister: Qubit[])
    : Unit is Adj + Ctl {
        PrepareQuantumROMState(nBitsPrecision, nCoeffs, nBitsIndices, keepCoeff, altIndex, signs, indexRegister, [signQubit], garbageRegister);
    }

    // Used in QuantumROM implementation.
    internal function QuantumROMBitStringWriterByIndex(idx : Int, keepCoeff : Int[], altIndex : Int[], data : Int[])
    : ((LittleEndian, LittleEndian, LittleEndian, LittleEndian) => Unit is Adj + Ctl) {
        return WriteQuantumROMBitString(idx, keepCoeff, altIndex, data, _, _, _, _);
    }

    // Used in QuantumROM implementation.
    internal operation WriteQuantumROMBitString(idx: Int, keepCoeff: Int[], altIndex: Int[], data : Int[], keepCoeffRegister: LittleEndian, altIndexRegister: LittleEndian, dataRegister : LittleEndian, altDataRegister : LittleEndian)
    : Unit is Adj + Ctl {
        ApplyXorInPlace(keepCoeff[idx], keepCoeffRegister);
        ApplyXorInPlace(altIndex[idx], altIndexRegister);
        if (Length(dataRegister!) > 0) {
            ApplyXorInPlace(data[idx], dataRegister);
            ApplyXorInPlace(data[altIndex[idx]], altDataRegister);
        }
    }

}
