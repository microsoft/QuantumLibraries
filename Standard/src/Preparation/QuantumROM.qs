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
    /// Returns an operation that prepares a a purification of a given mixed state. 
    /// A "purified mixed state" refers to states of the form |ψ⟩ = Σᵢ √𝑝ᵢ |𝑖⟩ |garbageᵢ⟩ specified by a vector of 
    /// coefficients {𝑝ᵢ}. States of this form can be reduced to mixed states ρ ≔ 𝑝ᵢ |𝑖⟩⟨𝑖| by tracing over the "garbage" 
    /// register (that is, a mixed state that is diagonal in the computational basis).
    /// 
    /// See https://arxiv.org/pdf/1805.03662.pdf?page=15 for further discussion.
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
    ///     \rho = \sum_{j = 0}^{N-1} \frac{|alpha_j|}{\sum_k |\alpha_k|} \ket{j}\bra{j},
    /// \end{align}
    /// $$
    /// where each $p_j$ is an approximation to the given coefficient $\alpha_j$
    /// such that
    /// $$
    /// \begin{align}
    ///     \left| p_j - \frac{ |\alpha_j| }{ \sum_k |\alpha_k| } \right| \le \frac{\epsilon}{N}
    /// \end{align}
    /// $$
    /// for each $j$.
    ///
    /// When passed an index register and a register of garbage qubits,
    /// initially in the state $\ket{0} \ket{00\cdots 0}$, the returned operation
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
    /// ```qsharp
    /// let coefficients = [1.0, 2.0, 3.0, 4.0, 5.0];
    /// let targetError = 1e-3;
    /// let purifiedState = PurifiedMixedState(targetError, coefficients);
    /// using (indexRegister = Qubit[purifiedState::Requirements::NIndexQubits]) {
    ///     using (garbageRegister = Qubit[purifiedState::Requirements::NGarbageQubits]) {
    ///         purifiedState::Prepare(LittleEndian(indexRegister), [], garbageRegister);
    ///     }
    /// }
    /// ```
    ///
    /// # References
    /// - Encoding Electronic Spectra in Quantum Circuits with Linear T Complexity
    ///   Ryan Babbush, Craig Gidney, Dominic W. Berry, Nathan Wiebe, Jarrod McClean, Alexandru Paler, Austin Fowler, Hartmut Neven
    ///   https://arxiv.org/abs/1805.03662
    ///
    /// # See Also
    /// - Microsoft.Quantum.Preparation.PurifiedMixedStateWithData
    function PurifiedMixedState(targetError : Double, coefficients : Double[])
    : MixedStatePreparation {
        let nBitsPrecision = -Ceiling(Lg(0.5 * targetError)) + 1;
        let positiveCoefficients = Mapped(AbsD, coefficients);
        let (oneNorm, keepCoeff, altIndex) = _QuantumROMDiscretization(nBitsPrecision, positiveCoefficients);
        let nCoeffs = Length(positiveCoefficients);
        let nBitsIndices = Ceiling(Lg(IntAsDouble(nCoeffs)));

        let op = PrepareQuantumROMState(nBitsPrecision, nCoeffs, nBitsIndices, keepCoeff, altIndex, [], _, _, _);
        let qubitCounts = PurifiedMixedStateRequirements(targetError, nCoeffs);
        return MixedStatePreparation(qubitCounts, oneNorm, op);
    }

    internal function SplitSign(coefficient : Double) : (Double, Int) {
        return (AbsD(coefficient), coefficient < 0.0 ? 1 | 0);
    }

    /// # Summary
    /// Returns an operation that prepares a a purification of a given mixed
    /// state, entangled with a register representing a given collection of data.
    /// A "purified mixed state with data" refers to a state of the form Σᵢ √𝑝ᵢ |𝑖⟩ |𝑥ᵢ⟩ |garbageᵢ⟩, 
    /// where each 𝑥ᵢ is a bitstring encoding additional data associated with the register |𝑖⟩.
    ///
    /// See https://arxiv.org/pdf/1805.03662.pdf?page=15 for further discussion.
    ///
    /// # Description
    /// Uses the Quantum ROM technique to represent a given density matrix,
    /// returning that representation as a state preparation operation.
    ///
    /// In particular, given a list of $N$ coefficients $\alpha_j$, and a
    /// bitstring $\vec{x}_j$ associated with each coefficient, this
    /// function returns an operation that uses the Quantum ROM technique to
    /// prepare an approximation
    /// $$
    /// \begin{align}
    ///     \tilde\rho = \sum_{j = 0}^{N - 1} p_j \ket{j}\bra{j} \otimes \ket{\vec{x}_j}\bra{\vec{x}_j}
    /// \end{align}
    /// $$
    /// of the mixed state
    /// $$
    /// \begin{align}
    ///     \rho = \sum_{j = 0}^{N-1}\ frac{|alpha_j|}{\sum_k |\alpha_k|} \ket{j}\bra{j} \otimes \ket{\vec{x}_j}\bra{\vec{x}_j},
    /// \end{align}
    /// $$
    /// where each $p_j$ is an approximation to the given coefficient $\alpha_j$
    /// such that
    /// $$
    /// \begin{align}
    ///     \left| p_j - \frac{ |\alpha_j| }{ \sum_k |\alpha_k| } \right| \le \frac{\epsilon}{N}
    /// \end{align}
    /// $$
    /// for each $j$.
    ///
    /// When passed an index register and a register of garbage qubits,
    /// initially in the state $\ket{0} \ket{00\cdots 0}$, the returned operation
    /// prepares both registers into the purification of $\tilde \rho$,
    /// $$
    /// \begin{align}
    ///     \sum_{j=0}^{N-1} \sqrt{p_j} \ket{j} \ket{\vec{x}_j} \ket{\text{garbage}_j},
    /// \end{align}
    /// $$
    /// such that resetting and deallocating the garbage register enacts the
    /// desired preparation to within the target error $\epsilon$.
    ///
    /// # Input
    /// ## targetError
    /// The target error $\epsilon$.
    /// ## coefficients
    /// Array of $N$ coefficients specifying the probability of basis states,
    /// along with the bitstring $\vec{x}_j$ associated with each coefficient.
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
    /// # References
    /// - Encoding Electronic Spectra in Quantum Circuits with Linear T Complexity
    ///   Ryan Babbush, Craig Gidney, Dominic W. Berry, Nathan Wiebe, Jarrod McClean, Alexandru Paler, Austin Fowler, Hartmut Neven
    ///   https://arxiv.org/abs/1805.03662
    ///
    /// # See Also
    /// - Microsoft.Quantum.Preparation.PurifiedMixedState
    function PurifiedMixedStateWithData(targetError : Double, coefficients : (Double, Bool[])[]) : MixedStatePreparation {
        let nBitsPrecision = -Ceiling(Lg(0.5 * targetError)) + 1;
        let positiveCoefficients = Mapped(AbsD, Mapped(Fst, coefficients));
        let data = Mapped(Snd, coefficients);
        let (oneNorm, keepCoeff, altIndex) = _QuantumROMDiscretization(nBitsPrecision, positiveCoefficients);
        let nCoeffs = Length(positiveCoefficients);
        let nBitsIndices = Ceiling(Lg(IntAsDouble(nCoeffs)));

        let op = PrepareQuantumROMState(nBitsPrecision, nCoeffs, nBitsIndices, keepCoeff, altIndex, data, _, _, _);
        let qubitCounts = PurifiedMixedStateRequirements(targetError, nCoeffs);
        return MixedStatePreparation(qubitCounts w/ NGarbageQubits <- qubitCounts::NGarbageQubits + 1, oneNorm, op);
    }

    /// # Summary
    /// Returns the total number of qubits that must be allocated
    /// in order to apply the operation returned by
    /// @"Microsoft.Quantum.Preparation.PurifiedMixedState".
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
    /// @"Microsoft.Quantum.Preparation.PurifiedMixedState" function.
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
            RoundedDiscretizationCoefficients(_, oneNorm, nCoefficients, barHeight),
            coefficients
        );

        // Calculate difference between number of discretized bars vs. maximum
        mutable bars = 0;
        for idxCoeff in IndexRange(keepCoeff) {
            set bars += keepCoeff[idxCoeff] - barHeight;
        }

        // Uniformly distribute excess bars across coefficients.
        for idx in 0..AbsI(bars) - 1 {
            set keepCoeff w/= idx <- keepCoeff[idx] + (bars > 0 ? -1 | +1);
        }

        mutable barSink = [];
        mutable barSource = [];

        for idxCoeff in IndexRange(keepCoeff) {
            if (keepCoeff[idxCoeff] > barHeight) {
                set barSource += [idxCoeff];
            } elif (keepCoeff[idxCoeff] < barHeight) {
                set barSink += [idxCoeff];
            }
        }

        for rep in 0..nCoefficients * 10 {
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
    internal function RoundedDiscretizationCoefficients(coefficient: Double, oneNorm: Double, nCoefficients: Int, barHeight: Int)
    : Int {
        return Round((AbsD(coefficient) / oneNorm) * IntAsDouble(nCoefficients) * IntAsDouble(barHeight));
    }

    // Used in QuantumROM implementation.
    internal operation PrepareQuantumROMState(
        nBitsPrecision: Int, nCoeffs: Int, nBitsIndices: Int,
        keepCoeff: Int[], altIndex: Int[], data : Bool[][],
        indexRegister: LittleEndian, dataQubits : Qubit[], garbageRegister: Qubit[]
    )
    : Unit is Adj + Ctl {
        let garbageIdx0 = nBitsIndices;
        let garbageIdx1 = garbageIdx0 + nBitsPrecision;
        let garbageIdx2 = garbageIdx1 + nBitsPrecision;
        let garbageIdx3 = garbageIdx2 + 1;

        let altIndexRegister = LittleEndian(garbageRegister[0..garbageIdx0 - 1]);
        let keepCoeffRegister = LittleEndian(garbageRegister[garbageIdx0..garbageIdx1 - 1]);
        let uniformKeepCoeffRegister = LittleEndian(garbageRegister[garbageIdx1..garbageIdx2 - 1]);
        let flagQubit = garbageRegister[garbageIdx3 - 1];
        let dataRegister = dataQubits;
        let altDataRegister = garbageRegister[garbageIdx3...];

        // Create uniform superposition over index and alt coeff register.
        PrepareUniformSuperposition(nCoeffs, indexRegister);
        ApplyToEachCA(H, uniformKeepCoeffRegister!);

        // Write bitstrings to altIndex and keepCoeff register.
        let unitaryGenerator = (nCoeffs, QuantumROMBitStringWriterByIndex(_, keepCoeff, altIndex, data));
        MultiplexOperationsFromGenerator(unitaryGenerator, indexRegister, (keepCoeffRegister, altIndexRegister, dataRegister, altDataRegister));

        // Perform comparison
        CompareUsingRippleCarry(uniformKeepCoeffRegister, keepCoeffRegister, flagQubit);

        let indexRegisterSize = Length(indexRegister!);

        // Swap in register based on comparison
        ApplyToEachCA((Controlled SWAP)([flagQubit], _), Zipped(indexRegister! + dataRegister, altIndexRegister! + altDataRegister));
    }

    // Used in QuantumROM implementation.
    internal function QuantumROMBitStringWriterByIndex(idx : Int, keepCoeff : Int[], altIndex : Int[], data : Bool[][])
    : ((LittleEndian, LittleEndian, Qubit[], Qubit[]) => Unit is Adj + Ctl) {
        return WriteQuantumROMBitString(idx, keepCoeff, altIndex, data, _, _, _, _);
    }

    // Used in QuantumROM implementation.
    internal operation WriteQuantumROMBitString(idx: Int, keepCoeff: Int[], altIndex: Int[], data : Bool[][], keepCoeffRegister: LittleEndian, altIndexRegister: LittleEndian, dataRegister : Qubit[], altDataRegister : Qubit[])
    : Unit is Adj + Ctl {
        if (keepCoeff[idx] >= 0) {
            ApplyXorInPlace(keepCoeff[idx], keepCoeffRegister);
        }
        ApplyXorInPlace(altIndex[idx], altIndexRegister);
        if (Length(dataRegister) > 0) {
            ApplyToEachCA(CControlledCA(X), Zipped(data[idx], dataRegister));
            ApplyToEachCA(CControlledCA(X), Zipped(data[altIndex[idx]], altDataRegister));
        }
    }

}
