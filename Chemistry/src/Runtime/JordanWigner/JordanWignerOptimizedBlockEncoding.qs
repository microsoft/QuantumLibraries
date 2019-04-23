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
    function _GetOptimizedBETermIndexCoeff_ (term : OptimizedBETermIndex) : Double {
        
        let (a, b, c, d, e, f) = term!;
        return a;
    }
    
    
    // Get OptimizedBEGeneratorSystem coefficients
    function _OptimizedBEGeneratorSystemCoeff_ (optimizedBEGeneratorSystem : OptimizedBEGeneratorSystem) : Double[] {
        
        let (nTerms, oneNorm, intToGenIdx) = optimizedBEGeneratorSystem!;
        mutable coefficients = new Double[nTerms];
        
        for (idx in 0 .. nTerms - 1) {
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
    function _ZTermToPauliMajIdx_ (term : GeneratorIndex) : OptimizedBETermIndex {
        
        let ((idxTermType, coeff), idxFermions) = term!;
        mutable signQubit = false;
        
        if (coeff[0] < 0.0) {
            set signQubit = true;
        }
        
        let selectZControlRegisters = [true];
        let OptimizedBEControlRegisters = new Bool[0];
        let pauliBases = new Int[0];
        let indexRegisters = idxFermions;
        return OptimizedBETermIndex(coeff[0], signQubit, selectZControlRegisters, OptimizedBEControlRegisters, pauliBases, indexRegisters);
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
    function _ZZTermToPauliMajIdx_ (term : GeneratorIndex) : OptimizedBETermIndex {
        
        let ((idxTermType, coeff), idxFermions) = term!;
        mutable signQubit = false;
        
        if (coeff[0] < 0.0) {
            set signQubit = true;
        }
        
        let selectZControlRegisters = [true, true];
        let OptimizedBEControlRegisters = new Bool[0];
        let pauliBases = new Int[0];
        let indexRegisters = idxFermions;
        return OptimizedBETermIndex(coeff[0], signQubit, selectZControlRegisters, OptimizedBEControlRegisters, pauliBases, indexRegisters);
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
    function _PQTermToPauliMajIdx_ (term : GeneratorIndex) : OptimizedBETermIndex {
        
        let ((idxTermType, coeff), idxFermions) = term!;
        mutable signQubit = false;
        
        if (coeff[0] < 0.0) {
            set signQubit = true;
        }
        
        let selectZControlRegisters = new Bool[0];
        let OptimizedBEControlRegisters = [true, true];
        let pauliBases = [1, 2];
        let indexRegisters = idxFermions;
        return OptimizedBETermIndex(2.0 * coeff[0], signQubit, selectZControlRegisters, OptimizedBEControlRegisters, pauliBases, indexRegisters);
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
    function _PQandPQQRTermToPauliMajIdx_ (term : GeneratorIndex) : OptimizedBETermIndex {
        
        let ((idxTermType, coeff), idxFermions) = term!;
        mutable signQubit = false;
        
        if (coeff[0] < 0.0) {
            set signQubit = true;
        }
        
        if (Length(idxFermions) == 2) {
            return _PQTermToPauliMajIdx_(term);
        }
        else {
            let qubitPidx = idxFermions[0];
            let qubitQidx = idxFermions[1];
            let qubitRidx = idxFermions[3];
            let selectZControlRegisters = [false, true];
            let OptimizedBEControlRegisters = [true, false, true];
            let pauliBases = [1, 2];
            let indexRegisters = [qubitPidx, qubitQidx, qubitRidx];
            return OptimizedBETermIndex(2.0 * coeff[0], signQubit, selectZControlRegisters, OptimizedBEControlRegisters, pauliBases, indexRegisters);
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
        let OptimizedBEControlRegisters = [true, true, true, true];
        let indexRegisters = idxFermions;
        
        for (idxOp in 0 .. 3) {
            
            if (IsNotZero(v0123[idxOp])) {
                mutable signQubit = false;
                
                if (v0123[idxOp] < 0.0) {
                    set signQubit = true;
                }
                
                let newCoeff = (2.0 * 0.25) * v0123[idxOp];
                set majIdxes w/= nonZero <- OptimizedBETermIndex(newCoeff, signQubit, selectZControlRegisters, OptimizedBEControlRegisters, ops[idxOp], indexRegisters);
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
        
        for (idx in IndexRange(ZData)) {
            
            // Array of Arrays of Length 1
            set majIdxes w/= idx <- _ZTermToPauliMajIdx_(HTermToGenIdx(ZData[idx], [0]));
        }
        
        set startIdx = Length(ZData);
        
        for (idx in IndexRange(ZZData)) {
            
            // Array of Arrays of Length 1
            set majIdxes w/= startIdx + idx <- _ZZTermToPauliMajIdx_(HTermToGenIdx(ZZData[idx], [1]));
        }
        
        set startIdx = startIdx + Length(ZZData);
        
        for (idx in IndexRange(PQandPQQRData)) {
            
            // Array of Arrays of Length 1
            set majIdxes w/= startIdx + idx <- _PQandPQQRTermToPauliMajIdx_(HTermToGenIdx(PQandPQQRData[idx], [2]));
        }
        
        set startIdx = startIdx + Length(PQandPQQRData);
        mutable finalIdx = startIdx;
        
        for (idx in 0 .. Length(h0123Data) - 1) {
            
            // Array of Arrays of Length up to 4
            let genArr = _V0123TermToPauliMajIdx_(HTermToGenIdx(h0123Data[idx], [3]));
            
            for (idx0123 in IndexRange(genArr)) {
                set majIdxes w/= finalIdx <- genArr[idx0123];
                set finalIdx = finalIdx + 1;
            }
        }
        
        mutable oneNorm = 0.0;
        
        for (idx in 0 .. finalIdx - 1) {
            set oneNorm = oneNorm + AbsD(_GetOptimizedBETermIndexCoeff_(majIdxes[idx]));
        }
        
        return OptimizedBEGeneratorSystem(finalIdx, oneNorm, LookupFunction(majIdxes[0 .. finalIdx - 1]));
    }
    
    
    operation _ToJordanWignerSelectInput__ (idx : Int, optimizedBEGeneratorSystem : OptimizedBEGeneratorSystem, signQubit : Qubit, selectZControlRegisters : Qubit[], OptimizedBEControlRegisters : Qubit[], pauliBasesIdx : LittleEndian, indexRegisters : LittleEndian[]) : Unit {
        
        body (...) {
            let (nTerms, oneNorm, intToGenIdx) = optimizedBEGeneratorSystem!;
            let (coeff, signQubitSet, selectZControlRegistersSet, OptimizedBEControlRegistersSet, pauliBasesSet, indexRegistersSet) = (intToGenIdx(idx))!;
            
            // Write bit to apply - signQubit
            if (signQubitSet == true) {
                X(signQubit);
            }
            
            // Write bit to activate selectZ operator
            for (idxSet in IndexRange(selectZControlRegistersSet)) {
                
                if (selectZControlRegistersSet[idxSet] == true) {
                    X(selectZControlRegisters[idxSet]);
                }
            }
            
            // Write bit to activate OptimizedBEXY operator
            for (idxSet in IndexRange(OptimizedBEControlRegistersSet)) {
                
                if (OptimizedBEControlRegistersSet[idxSet] == true) {
                    X(OptimizedBEControlRegisters[idxSet]);
                }
            }
            
            // Write bitstring to apply desired XZ... or YZ... Pauli string
            for (idxSet in IndexRange(indexRegistersSet)) {
                ApplyXorInPlace(indexRegistersSet[idxSet], indexRegisters[idxSet]);
            }
            
            // Crete state to select uniform superposition of X and Y operators.
            if (Length(pauliBasesSet) == 2) {
                
                // for PQ or PQQR terms, create |00> + |11>
                ApplyXorInPlace(0, pauliBasesIdx);
            }
            elif (Length(pauliBasesSet) == 4) {
                
                // for PQRS terms, create |abcd> + |a^ b^ c^ d^>
                if (pauliBasesSet[2] == 1 and pauliBasesSet[3] == 1) {
                    ApplyXorInPlace(1, pauliBasesIdx);
                }
                elif (pauliBasesSet[2] == 2 and pauliBasesSet[3] == 2) {
                    ApplyXorInPlace(2, pauliBasesIdx);
                }
                elif (pauliBasesSet[2] == 1 and pauliBasesSet[3] == 2) {
                    ApplyXorInPlace(3, pauliBasesIdx);
                }
                elif (pauliBasesSet[2] == 2 and pauliBasesSet[3] == 1) {
                    ApplyXorInPlace(4, pauliBasesIdx);
                }
            }
        }
        
        adjoint invert;
        controlled distribute;
        controlled adjoint distribute;
    }
    
    
    function _ToJordanWignerSelectInput_ (idx : Int, optimizedBEGeneratorSystem : OptimizedBEGeneratorSystem) : ((Qubit, Qubit[], Qubit[], LittleEndian, LittleEndian[]) => Unit is Adj + Ctl) {
        
        return _ToJordanWignerSelectInput__(idx, optimizedBEGeneratorSystem, _, _, _, _, _);
    }
    
    
    operation _ToPauliBases__ (idx : Int, pauliBases : Qubit[]) : Unit {
        
        body (...) {
            let pauliBasesSet = [[1, 1, 1, 1], [1, 1, 2, 2], [1, 2, 1, 2], [1, 2, 2, 1]];
            H(pauliBases[0]);
            
            if (idx > 0) {
                
                for (idxSet in 1 .. Length(pauliBasesSet[0]) - 1) {
                    
                    if ((pauliBasesSet[idx - 1])[idxSet] == 2) {
                        X(pauliBases[idxSet]);
                    }
                    
                    CNOT(pauliBases[0], pauliBases[idxSet]);
                }
            }
        }
        
        adjoint invert;
        controlled distribute;
        controlled adjoint distribute;
    }
    
    
    function _ToPauliBases_ (idx : Int) : (Qubit[] => Unit is Adj + Ctl) {
        
        return _ToPauliBases__(idx, _);
    }
    
    
    // This prepares the state that selects _JordanWignerSelect_;
    operation _JordanWignerOptimizedBlockEncodingStatePrep__ (targetError : Double, optimizedBEGeneratorSystem : OptimizedBEGeneratorSystem, qROMIdxRegister : LittleEndian, qROMGarbage : Qubit[], signQubit : Qubit, selectZControlRegisters : Qubit[], OptimizedBEControlRegisters : Qubit[], pauliBases : Qubit[], pauliBasesIdx : LittleEndian, indexRegisters : LittleEndian[]) : Unit {
        
        body (...) {
            let (nTerms, oneNorm0, intToGenIdx) = optimizedBEGeneratorSystem!;
            let coefficients = _OptimizedBEGeneratorSystemCoeff_(optimizedBEGeneratorSystem);
            let (qROMQubitCount, oneNorm, qROMUnitary) = QuantumROM(targetError, coefficients);
            let unitaryGenerator = (nTerms, _ToJordanWignerSelectInput_(_, optimizedBEGeneratorSystem));
            let pauliBasesUnitaryGenerator = (5, _ToPauliBases_);
            
            //let multiplexer = MultiplexerFromGenerator;
            qROMUnitary(qROMIdxRegister, qROMGarbage);
            MultiplexOperationsFromGenerator(unitaryGenerator, qROMIdxRegister, (signQubit, selectZControlRegisters, OptimizedBEControlRegisters, pauliBasesIdx, indexRegisters));
            MultiplexOperationsFromGenerator(pauliBasesUnitaryGenerator, pauliBasesIdx, pauliBases);
        }
        
        adjoint invert;
        controlled distribute;
        controlled adjoint distribute;
    }
    
    
    function _JordanWignerOptimizedBlockEncodingQubitManager_ (targetError : Double, nCoeffs : Int, nZ : Int, nMaj : Int, nIdxRegQubits : Int, ctrlRegister : Qubit[]) : ((LittleEndian, Qubit[], Qubit, Qubit[], Qubit[], Qubit[], LittleEndian, LittleEndian[]), (Qubit, Qubit[], Qubit[], Qubit[], LittleEndian[]), Qubit[]) {
        
        let ((qROMIdx, qROMGarbage), rest0) = _QuantumROMQubitManager(targetError, nCoeffs, ctrlRegister);
        let ((signQubit, selectZControlRegisters, optimizedBEControlRegisters, pauliBases, indexRegisters, tmp), rest1) = _JordanWignerSelectQubitManager_(nZ, nMaj, nIdxRegQubits, rest0, new Qubit[0]);
        let registers = Partitioned([3], rest1);
        let pauliBasesIdx = LittleEndian(registers[0]);
        return ((qROMIdx, qROMGarbage, signQubit, selectZControlRegisters, optimizedBEControlRegisters, pauliBases, pauliBasesIdx, indexRegisters), (signQubit, selectZControlRegisters, optimizedBEControlRegisters, pauliBases, indexRegisters), registers[1]);
    }
    
    function _JordanWignerOptimizedBlockEncodingQubitCount_ (targetError : Double, nCoeffs : Int, nZ : Int, nMaj : Int, nIdxRegQubits : Int, nTarget : Int) : ((Int, Int), (Int, Int, Int, Int, Int, Int, Int, Int[], Int)) {
        
        let (nSelectTotal, (a0, a1, a2, a3, a4)) = _JordanWignerSelectQubitCount_(nZ, nMaj, nIdxRegQubits);
        let (nQROMTotal, (b0, b1)) = QuantumROMQubitCount(targetError, nCoeffs);
        let pauliBasesIdx = 3;
        return (((nSelectTotal + nQROMTotal) + pauliBasesIdx, nTarget), (b0, b1, a0, a1, a2, a3, pauliBasesIdx, a4, nTarget));
    }
    
    
    operation _JordanWignerOptimizedBlockEncodingStatePrep_ (targetError : Double, nCoeffs : Int, optimizedBEGeneratorSystem : OptimizedBEGeneratorSystem, nZ : Int, nMaj : Int, nIdxRegQubits : Int, ctrlRegister : Qubit[]) : Unit {
        
        body (...) {
            let (statePrepRegister, selectRegister, rest) = _JordanWignerOptimizedBlockEncodingQubitManager_(targetError, nCoeffs, nZ, nMaj, nIdxRegQubits, ctrlRegister);
            let statePrepOp = _JordanWignerOptimizedBlockEncodingStatePrep__(targetError, optimizedBEGeneratorSystem, _, _, _, _, _, _, _, _);
            statePrepOp(statePrepRegister);
        }
        
        adjoint invert;
        controlled distribute;
        controlled adjoint distribute;
    }
    
    
    operation _JordanWignerOptimizedBlockEncodingSelect_ (targetError : Double, nCoeffs : Int, optimizedBEGeneratorSystem : OptimizedBEGeneratorSystem, nZ : Int, nMaj : Int, nIdxRegQubits : Int, ctrlRegister : Qubit[], targetRegister : Qubit[]) : Unit {
        
        body (...) {
            let (statePrepRegister, selectRegister, rest) = _JordanWignerOptimizedBlockEncodingQubitManager_(targetError, nCoeffs, nZ, nMaj, nIdxRegQubits, ctrlRegister);
            let selectOp = _JordanWignerSelect_(_, _, _, _, _, targetRegister);
            selectOp(selectRegister);
        }
        
        adjoint invert;
        controlled distribute;
        controlled adjoint distribute;
    }
    
    
    function _JordanWignerOptimizedBlockEncoding_ (targetError : Double, data : JWOptimizedHTerms, nSpinOrbitals : Int) : ((Int, Int), (Double, BlockEncodingReflection)) {
        
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
    
    
    function _JordanWignerOptimizedQuantumWalkByQubitization_ (targetError : Double, data : JWOptimizedHTerms, nSpinOrbitals : Int) : ((Int, Int), (Double, ((Qubit[], Qubit[]) => Unit is Adj + Ctl))) {
        
        let ((nCtrlRegisterQubits, nTargetRegisterQubits), (oneNorm, blockEncodingReflection)) = _JordanWignerOptimizedBlockEncoding_(targetError, data, nSpinOrbitals);
        return ((nCtrlRegisterQubits, nTargetRegisterQubits), (oneNorm, QuantumWalkByQubitization(blockEncodingReflection)));
    }
    
}


