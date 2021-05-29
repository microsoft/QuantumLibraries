// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Chemistry.JordanWigner {
    open Microsoft.Quantum.Simulation;
    open Microsoft.Quantum.Arithmetic;
    open Microsoft.Quantum.Preparation;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Convert;
    open Microsoft.Quantum.Chemistry;
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Math;

    /// # Summary
    /// Term data in the optimized block-encoding algorithm.
    newtype OptimizedBETermIndex = (Double, Bool, Bool[], Bool[], Int[], Int[]);

    /// # Summary
    /// Function that returns `OptimizedBETermIndex` data for term `n` given an
    /// integer `n`, together with the number of terms in the first `Int` and
    /// the sum of absolute-values of all term coefficients in the `Double`.
    newtype OptimizedBEGeneratorSystem = (Int, Double, (Int -> OptimizedBETermIndex));


    // Get OptimizedBETermIndex coefficients
    internal function _GetOptimizedBETermIndexCoeff_ (term : OptimizedBETermIndex) : Double {
        let (a, b, c, d, e, f) = term!;
        return a;
    }


    // Get OptimizedBEGeneratorSystem coefficients
    internal function _OptimizedBEGeneratorSystemCoeff_ (optimizedBEGeneratorSystem : OptimizedBEGeneratorSystem) : Double[] {

        let (nTerms, oneNorm, intToGenIdx) = optimizedBEGeneratorSystem!;
        mutable coefficients = new Double[nTerms];

        for idx in 0 .. nTerms - 1 {
            set coefficients w/= idx <- _GetOptimizedBETermIndexCoeff_(intToGenIdx(idx));
        }

        return coefficients;
    }


    /// # Summary
    /// Converts a GeneratorIndex describing a Z term to
    /// an expression 'GeneratorIndex[]' in terms of Paulis.
    ///
    /// # Input
    /// ## term
    /// `GeneratorIndex` representing a Z term.
    ///
    /// # Output
    /// 'GeneratorIndex[]' expressing Z term as Pauli terms.
    internal function _ZTermToPauliMajIdx_ (term : GeneratorIndex) : OptimizedBETermIndex {

        let ((idxTermType, coeff), idxFermions) = term!;
        mutable signQubit = false;

        if (coeff[0] < 0.0) {
            set signQubit = true;
        }

        let selectZControlRegisters = [true];
        let optimizedBEControlRegisters = new Bool[0];
        let pauliBases = new Int[0];
        let indexRegisters = idxFermions;
        return OptimizedBETermIndex(coeff[0], signQubit, selectZControlRegisters, optimizedBEControlRegisters, pauliBases, indexRegisters);
    }


    /// # Summary
    /// Converts a GeneratorIndex describing a ZZ term to
    /// an expression 'GeneratorIndex[]' in terms of Paulis.
    ///
    /// # Input
    /// ## term
    /// `GeneratorIndex` representing a ZZ term.
    ///
    /// # Output
    /// 'GeneratorIndex[]' expressing ZZ term as Pauli terms.
    internal function _ZZTermToPauliMajIdx_(term : GeneratorIndex) : OptimizedBETermIndex {

        let ((idxTermType, coeff), idxFermions) = term!;
        mutable signQubit = false;

        if coeff[0] < 0.0 {
            set signQubit = true;
        }

        let selectZControlRegisters = [true, true];
        let optimizedBEControlRegisters = new Bool[0];
        let pauliBases = new Int[0];
        let indexRegisters = idxFermions;
        return OptimizedBETermIndex(coeff[0], signQubit, selectZControlRegisters, optimizedBEControlRegisters, pauliBases, indexRegisters);
    }


    /// # Summary
    /// Converts a GeneratorIndex describing a PQ term to
    /// an expression 'GeneratorIndex[]' in terms of Paulis
    ///
    /// # Input
    /// ## term
    /// `GeneratorIndex` representing a PQ term.
    ///
    /// # Output
    /// 'GeneratorIndex[]' expressing PQ term as Pauli terms.
    internal function _PQTermToPauliMajIdx_ (term : GeneratorIndex) : OptimizedBETermIndex {

        let ((idxTermType, coeff), idxFermions) = term!;
        let sign = coeff[0] < 0.0;

        let selectZControlRegisters = new Bool[0];
        let optimizedBEControlRegisters = [true, true];
        let pauliBases = [1, 2];
        let indexRegisters = idxFermions;
        return OptimizedBETermIndex(2.0 * coeff[0], sign, selectZControlRegisters, optimizedBEControlRegisters, pauliBases, indexRegisters);
    }


    /// # Summary
    /// Converts a GeneratorIndex describing a PQ or PQQR term to
    /// an expression 'GeneratorIndex[]' in terms of Paulis
    ///
    /// # Input
    /// ## term
    /// `GeneratorIndex` representing a PQ or PQQR term.
    ///
    /// # Output
    /// 'GeneratorIndex[]' expressing PQ or PQQR term as Pauli terms.
    internal function _PQandPQQRTermToPauliMajIdx_ (term : GeneratorIndex) : OptimizedBETermIndex {
        let ((idxTermType, coeff), idxFermions) = term!;
        let sign = coeff[0] < 0.0;

        if Length(idxFermions) == 2 {
            return _PQTermToPauliMajIdx_(term);
        } else {
            let qubitPidx = idxFermions[0];
            let qubitQidx = idxFermions[1];
            let qubitRidx = idxFermions[3];
            let selectZControlRegisters = [false, true];
            let optimizedBEControlRegisters = [true, false, true];
            let pauliBases = [1, 2];
            let indexRegisters = [qubitPidx, qubitQidx, qubitRidx];
            return OptimizedBETermIndex(2.0 * coeff[0], sign, selectZControlRegisters, optimizedBEControlRegisters, pauliBases, indexRegisters);
        }
    }


    /// # Summary
    /// Converts a GeneratorIndex describing a PQRS term to
    /// an expression 'GeneratorIndex[]' in terms of Paulis
    ///
    /// # Input
    /// ## term
    /// `GeneratorIndex` representing a PQRS term.
    ///
    /// # Output
    /// 'GeneratorIndex[]' expressing PQRS term as Pauli terms.
    function _V0123TermToPauliMajIdx_ (term : GeneratorIndex) : OptimizedBETermIndex[] {

        let ((idxTermType, v0123), idxFermions) = term!;
        let qubitsPQ = idxFermions[0 .. 1];
        let qubitsRS = idxFermions[2 .. 3];
        let qubitsPQJW = RangeAsIntArray(qubitsPQ[0] + 1 .. qubitsPQ[1] - 1);
        let qubitsRSJW = RangeAsIntArray(qubitsRS[0] + 1 .. qubitsRS[1] - 1);
        let ops = [[1, 1, 1, 1], [1, 1, 2, 2], [1, 2, 1, 2], [1, 2, 2, 1], [2, 2, 2, 2], [2, 2, 1, 1], [2, 1, 2, 1], [2, 1, 1, 2]];
        mutable majIdxes = new OptimizedBETermIndex[4];
        mutable nonZero = 0;
        let selectZControlRegisters = new Bool[0];
        let optimizedBEControlRegisters = [true, true, true, true];
        let indexRegisters = idxFermions;

        for idxOp in 0 .. 3 {
            if IsNotZero(v0123[idxOp]) {
                let newCoeff = (2.0 * 0.25) * v0123[idxOp];
                set majIdxes w/= nonZero <- OptimizedBETermIndex(newCoeff, v0123[idxOp] < 0.0, selectZControlRegisters, optimizedBEControlRegisters, ops[idxOp], indexRegisters);
                set nonZero = nonZero + 1;
            }
        }

        return majIdxes[0 .. nonZero - 1];
    }


    /// # Summary
    /// Converts a Hamiltonian described by `JWOptimizedHTerms`
    /// to a `GeneratorSystem` expressed in terms of the Pauli
    /// `GeneratorIndex`.
    ///
    /// # Input
    /// ## data
    /// Description of Hamiltonian in `JWOptimizedHTerms` format.
    ///
    /// # Output
    /// Representation of Hamiltonian as `GeneratorSystem`.
    function OptimizedBlockEncodingGeneratorSystem (data : JWOptimizedHTerms) : OptimizedBEGeneratorSystem {

        let (ZData, ZZData, PQandPQQRData, h0123Data) = data!;
        mutable majIdxes = new OptimizedBETermIndex[((Length(ZData) + Length(ZZData)) + Length(PQandPQQRData)) + 4 * Length(h0123Data)];
        mutable startIdx = 0;

        for idx in IndexRange(ZData) {
            // Array of Arrays of Length 1
            set majIdxes w/= idx <- _ZTermToPauliMajIdx_(HTermToGenIdx(ZData[idx], [0]));
        }

        set startIdx = Length(ZData);

        for idx in IndexRange(ZZData) {
            // Array of Arrays of Length 1
            set majIdxes w/= startIdx + idx <- _ZZTermToPauliMajIdx_(HTermToGenIdx(ZZData[idx], [1]));
        }

        set startIdx = startIdx + Length(ZZData);

        for idx in IndexRange(PQandPQQRData) {

            // Array of Arrays of Length 1
            set majIdxes w/= startIdx + idx <- _PQandPQQRTermToPauliMajIdx_(HTermToGenIdx(PQandPQQRData[idx], [2]));
        }

        set startIdx = startIdx + Length(PQandPQQRData);
        mutable finalIdx = startIdx;

        for idx in 0 .. Length(h0123Data) - 1 {

            // Array of Arrays of Length up to 4
            let genArr = _V0123TermToPauliMajIdx_(HTermToGenIdx(h0123Data[idx], [3]));

            for idx0123 in IndexRange(genArr) {
                set majIdxes w/= finalIdx <- genArr[idx0123];
                set finalIdx = finalIdx + 1;
            }
        }

        mutable oneNorm = 0.0;

        for idx in 0 .. finalIdx - 1 {
            set oneNorm = oneNorm + AbsD(_GetOptimizedBETermIndexCoeff_(majIdxes[idx]));
        }

        return OptimizedBEGeneratorSystem(finalIdx, oneNorm, LookupFunction(majIdxes[0 .. finalIdx - 1]));
    }


    internal operation _ToJordanWignerSelectInput (
        idx : Int,
        optimizedBEGeneratorSystem : OptimizedBEGeneratorSystem,
        signQubit : Qubit,
        selectZControlRegisters : Qubit[],
        optimizedBEControlRegisters : Qubit[],
        pauliBasesIdx : LittleEndian,
        indexRegisters : LittleEndian[]
    ) : Unit is Adj + Ctl {
        let (nTerms, oneNorm, intToGenIdx) = optimizedBEGeneratorSystem!;
        let (coeff, signQubitSet, selectZControlRegistersSet, OptimizedBEControlRegistersSet, pauliBasesSet, indexRegistersSet) = (intToGenIdx(idx))!;

        // Write bit to apply - signQubit
        if (signQubitSet == true) {
            X(signQubit);
        }

        // Write bit to activate selectZ operator
        ApplyToEachCA(CControlledCA(X), Zipped(selectZControlRegistersSet, selectZControlRegisters));

        // Write bit to activate OptimizedBEXY operator
        ApplyToEachCA(CControlledCA(X), Zipped(OptimizedBEControlRegistersSet, optimizedBEControlRegisters));

        // Write bitstring to apply desired XZ... or YZ... Pauli string
        ApplyToEachCA(ApplyXorInPlace, Zipped(indexRegistersSet, indexRegisters));

        // Crete state to select uniform superposition of X and Y operators.
        if Length(pauliBasesSet) == 2 {
            // for PQ or PQQR terms, create |00> + |11>
            ApplyXorInPlace(0, pauliBasesIdx);
        } elif Length(pauliBasesSet) == 4 {
            // for PQRS terms, create |abcd> + |a^ b^ c^ d^>
            if pauliBasesSet[2] == 1 and pauliBasesSet[3] == 1 {
                ApplyXorInPlace(1, pauliBasesIdx);
            }
            elif pauliBasesSet[2] == 2 and pauliBasesSet[3] == 2 {
                ApplyXorInPlace(2, pauliBasesIdx);
            }
            elif pauliBasesSet[2] == 1 and pauliBasesSet[3] == 2 {
                ApplyXorInPlace(3, pauliBasesIdx);
            }
            elif pauliBasesSet[2] == 2 and pauliBasesSet[3] == 1 {
                ApplyXorInPlace(4, pauliBasesIdx);
            }
        }
    }


    internal function _ToJordanWignerSelectInput_ (idx : Int, optimizedBEGeneratorSystem : OptimizedBEGeneratorSystem)
    : ((Qubit, Qubit[], Qubit[], LittleEndian, LittleEndian[]) => Unit is Adj + Ctl) {
        return _ToJordanWignerSelectInput(idx, optimizedBEGeneratorSystem, _, _, _, _, _);
    }


    internal operation _ToPauliBases(idx : Int, pauliBases : Qubit[]) : Unit is Adj + Ctl {
        let pauliBasesSet = [[1, 1, 1, 1], [1, 1, 2, 2], [1, 2, 1, 2], [1, 2, 2, 1]];
        H(pauliBases[0]);

        if idx > 0 {
            for idxSet in 1 .. Length(pauliBasesSet[0]) - 1 {
                if (pauliBasesSet[idx - 1])[idxSet] == 2 {
                    X(pauliBases[idxSet]);
                }

                CNOT(pauliBases[0], pauliBases[idxSet]);
            }
        }
    }

    // This prepares the state that selects _JordanWignerSelect_;
    internal operation _JordanWignerOptimizedBlockEncodingStatePrep(targetError : Double, optimizedBEGeneratorSystem : OptimizedBEGeneratorSystem, qROMIdxRegister : LittleEndian, qROMGarbage : Qubit[], signQubit : Qubit, selectZControlRegisters : Qubit[], optimizedBEControlRegisters : Qubit[], pauliBases : Qubit[], pauliBasesIdx : LittleEndian, indexRegisters : LittleEndian[])
    : Unit is Adj + Ctl {
        let (nTerms, _, _) = optimizedBEGeneratorSystem!;
        let coefficients = _OptimizedBEGeneratorSystemCoeff_(optimizedBEGeneratorSystem);
        let purifiedState = PurifiedMixedState(targetError, coefficients);
        let unitaryGenerator = (nTerms, _ToJordanWignerSelectInput_(_, optimizedBEGeneratorSystem));
        let pauliBasesUnitaryGenerator = (5, CurriedOpCA(_ToPauliBases));

        purifiedState::Prepare(qROMIdxRegister, [], qROMGarbage);
        MultiplexOperationsFromGenerator(unitaryGenerator, qROMIdxRegister, (signQubit, selectZControlRegisters, optimizedBEControlRegisters, pauliBasesIdx, indexRegisters));
        MultiplexOperationsFromGenerator(pauliBasesUnitaryGenerator, pauliBasesIdx, pauliBases);
    }

    internal function _JordanWignerOptimizedBlockEncodingQubitManager_ (targetError : Double, nCoeffs : Int, nZ : Int, nMaj : Int, nIdxRegQubits : Int, ctrlRegister : Qubit[]) : ((LittleEndian, Qubit[], Qubit, Qubit[], Qubit[], Qubit[], LittleEndian, LittleEndian[]), (Qubit, Qubit[], Qubit[], Qubit[], LittleEndian[]), Qubit[]) {
        let requirements = PurifiedMixedStateRequirements(targetError, nCoeffs);
        let parts = Partitioned([requirements::NIndexQubits, requirements::NGarbageQubits], ctrlRegister);
        let ((qROMIdx, qROMGarbage), rest0) = ((LittleEndian(parts[0]), parts[1]), parts[2]);
        let ((signQubit, selectZControlRegisters, optimizedBEControlRegisters, pauliBases, indexRegisters, tmp), rest1) = _JordanWignerSelectQubitManager_(nZ, nMaj, nIdxRegQubits, rest0, []);
        let registers = Partitioned([3], rest1);
        let pauliBasesIdx = LittleEndian(registers[0]);
        return ((qROMIdx, qROMGarbage, signQubit, selectZControlRegisters, optimizedBEControlRegisters, pauliBases, pauliBasesIdx, indexRegisters), (signQubit, selectZControlRegisters, optimizedBEControlRegisters, pauliBases, indexRegisters), registers[1]);
    }

    internal function _JordanWignerOptimizedBlockEncodingQubitCount_ (targetError : Double, nCoeffs : Int, nZ : Int, nMaj : Int, nIdxRegQubits : Int, nTarget : Int) : ((Int, Int), (Int, Int, Int, Int, Int, Int, Int, Int[], Int)) {

        let (nSelectTotal, (a0, a1, a2, a3, a4)) = _JordanWignerSelectQubitCount_(nZ, nMaj, nIdxRegQubits);
        let (nQROMTotal, (b0, b1)) = (PurifiedMixedStateRequirements(targetError, nCoeffs))!;
        let pauliBasesIdx = 3;
        return (((nSelectTotal + nQROMTotal) + pauliBasesIdx, nTarget), (b0, b1, a0, a1, a2, a3, pauliBasesIdx, a4, nTarget));
    }


    internal operation _JordanWignerOptimizedBlockEncodingStatePrep_ (targetError : Double, nCoeffs : Int, optimizedBEGeneratorSystem : OptimizedBEGeneratorSystem, nZ : Int, nMaj : Int, nIdxRegQubits : Int, ctrlRegister : Qubit[])
    : Unit is Adj + Ctl {
        let (statePrepRegister, selectRegister, rest) = _JordanWignerOptimizedBlockEncodingQubitManager_(targetError, nCoeffs, nZ, nMaj, nIdxRegQubits, ctrlRegister);
        let statePrepOp = _JordanWignerOptimizedBlockEncodingStatePrep(targetError, optimizedBEGeneratorSystem, _, _, _, _, _, _, _, _);
        statePrepOp(statePrepRegister);
    }


    internal operation _JordanWignerOptimizedBlockEncodingSelect_(targetError : Double, nCoeffs : Int, optimizedBEGeneratorSystem : OptimizedBEGeneratorSystem, nZ : Int, nMaj : Int, nIdxRegQubits : Int, ctrlRegister : Qubit[], targetRegister : Qubit[])
    : Unit is Adj + Ctl {
        let (statePrepRegister, selectRegister, rest) = _JordanWignerOptimizedBlockEncodingQubitManager_(targetError, nCoeffs, nZ, nMaj, nIdxRegQubits, ctrlRegister);
        let selectOp = _JordanWignerSelect_(_, _, _, _, _, targetRegister);
        selectOp(selectRegister);
    }


    function _JordanWignerOptimizedBlockEncoding_(targetError : Double, data : JWOptimizedHTerms, nSpinOrbitals : Int) : ((Int, Int), (Double, BlockEncodingReflection)) {
        let nZ = 2;
        let nMaj = 4;
        let optimizedBEGeneratorSystem = OptimizedBlockEncodingGeneratorSystem(data);
        let (nCoeffs, oneNorm, tmp) = optimizedBEGeneratorSystem!;
        let nIdxRegQubits = Ceiling(Lg(IntAsDouble(nSpinOrbitals)));
        let ((nCtrlRegisterQubits, nTargetRegisterQubits), rest) = _JordanWignerOptimizedBlockEncodingQubitCount_(targetError, nCoeffs, nZ, nMaj, nIdxRegQubits, nSpinOrbitals);
        let statePrepOp = _JordanWignerOptimizedBlockEncodingStatePrep_(targetError, nCoeffs, optimizedBEGeneratorSystem, nZ, nMaj, nIdxRegQubits, _);
        let selectOp = _JordanWignerOptimizedBlockEncodingSelect_(targetError, nCoeffs, optimizedBEGeneratorSystem, nZ, nMaj, nIdxRegQubits, _, _);
        let blockEncodingReflection = BlockEncodingReflection(BlockEncoding(BlockEncodingByLCU(statePrepOp, selectOp)));
        return ((nCtrlRegisterQubits, nTargetRegisterQubits), (oneNorm, blockEncodingReflection));
    }

    internal function _JordanWignerOptimizedQuantumWalkByQubitization_(targetError : Double, data : JWOptimizedHTerms, nSpinOrbitals : Int) : ((Int, Int), (Double, ((Qubit[], Qubit[]) => Unit is Adj + Ctl))) {
        let ((nCtrlRegisterQubits, nTargetRegisterQubits), (oneNorm, blockEncodingReflection)) = _JordanWignerOptimizedBlockEncoding_(targetError, data, nSpinOrbitals);
        return ((nCtrlRegisterQubits, nTargetRegisterQubits), (oneNorm, QuantumWalkByQubitization(blockEncodingReflection)));
    }

}


