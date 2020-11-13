// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Chemistry.JordanWigner {
    open Microsoft.Quantum.Simulation;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Chemistry;
    open Microsoft.Quantum.Arrays;

    // This evolution set runs off data optimized for a Jordan–Wigner encoding.
    // This collects terms Z, ZZ, PQandPQQR, hpqrs separately.
    // This only apples the needed hpqrs XXXX XXYY terms.
    // Operations here are expressed in terms of Exp([...])

    // Convention for GeneratorIndex = ((Int[],Double[]), Int[])
    // We index single Paulis as 0 for I, 1 for X, 2 for Y, 3 for Z.
    // We index Pauli strings with arrays of integers e.g. a = [3,1,1,2] for ZXXY.
    // We assume the zeroth element of Double[] is the angle of rotation
    // We index the qubits that Pauli strings act on with arrays of integers e.g. q = [2,4,5,8] for Z_2 X_4 X_5, Y_8
    // An example of a Pauli string GeneratorIndex is thus ((a,b), q)

    // Consider the Hamiltonian H = 0.1 XI + 0.2 IX + 0.3 ZY
    // Its GeneratorTerms are (([1],b),[0]), 0.1),  (([1],b),[1]), 0.2),  (([3,2],b),[0,1]), 0.3).

    /// # Summary
    /// Applies time-evolution by a Z term described by a `GeneratorIndex`.
    ///
    /// # Input
    /// ## term
    /// `GeneratorIndex` representing a Z term.
    /// ## stepSize
    /// Duration of time-evolution.
    /// ## qubits
    /// Qubits of Hamiltonian.
    operation _ApplyJordanWignerZTerm_ (term : GeneratorIndex, stepSize : Double, qubits : Qubit[]) : Unit {
        body (...) {
            let ((idxTermType, coeff), idxFermions) = term!;
            let angle = (1.0 * coeff[0]) * stepSize;
            let qubit = qubits[idxFermions[0]];
            Exp([PauliZ], angle, [qubit]);
        }
        
        adjoint invert;
        controlled distribute;
        controlled adjoint distribute;
    }
    
    
    /// # Summary
    /// Applies time-evolution by a ZZ term described by a `GeneratorIndex`.
    ///
    /// # Input
    /// ## term
    /// `GeneratorIndex` representing a ZZ term.
    /// ## stepSize
    /// Duration of time-evolution.
    /// ## qubits
    /// Qubits of Hamiltonian.
    operation _ApplyJordanWignerZZTerm_ (term : GeneratorIndex, stepSize : Double, qubits : Qubit[]) : Unit {
        
        body (...) {
            let ((idxTermType, coeff), idxFermions) = term!;
            let angle = (1.0 * coeff[0]) * stepSize;
            let qubitsZZ = Subarray(idxFermions[0 .. 1], qubits);
            Exp([PauliZ, PauliZ], angle, qubitsZZ);
        }
        
        adjoint invert;
        controlled distribute;
        controlled adjoint distribute;
    }
    
    
    /// # Summary
    /// Applies time-evolution by a PQ term described by a `GeneratorIndex`.
    ///
    /// # Input
    /// ## term
    /// `GeneratorIndex` representing a PQ term.
    /// ## stepSize
    /// Duration of time-evolution.
    /// ## extraParityQubits
    /// Optional parity qubits that flip the sign of time-evolution.
    /// ## qubits
    /// Qubits of Hamiltonian.
    operation _ApplyJordanWignerPQTerm_ (term : GeneratorIndex, stepSize : Double, extraParityQubits : Qubit[], qubits : Qubit[]) : Unit {
        
        body (...) {
            let ((idxTermType, coeff), idxFermions) = term!;
            let angle = (1.0 * coeff[0]) * stepSize;
            let qubitsPQ = Subarray(idxFermions[0 .. 1], qubits);
            let qubitsJW = qubits[idxFermions[0] + 1 .. idxFermions[1] - 1];
            let ops = [[PauliX, PauliX], [PauliY, PauliY]];
            
            for (idxOp in IndexRange(ops)) {
                Exp(ops[idxOp] + ConstantArray(Length(qubitsJW) + Length(extraParityQubits), PauliZ), angle, (qubitsPQ + qubitsJW) + extraParityQubits);
            }
        }
        
        adjoint invert;
        controlled distribute;
        controlled adjoint distribute;
    }
    
    
    /// # Summary
    /// Applies time-evolution by a PQ or PQQR term described by a `GeneratorIndex`.
    ///
    /// # Input
    /// ## term
    /// `GeneratorIndex` representing a PQ or PQQR term.
    /// ## stepSize
    /// Duration of time-evolution.
    /// ## qubits
    /// Qubits of Hamiltonian.
    operation _ApplyJordanWignerPQandPQQRTerm_ (term : GeneratorIndex, stepSize : Double, qubits : Qubit[]) : Unit {
        body (...) {
            let ((idxTermType, coeff), idxFermions) = term!;
            let angle = (1.0 * coeff[0]) * stepSize;
            let qubitQidx = idxFermions[1];
            
            // For all cases, do the same thing:
            // p < r < q (1/4)(1-Z_q)(Z_{r-1,p+1})(X_p X_r + Y_p Y_r) (same as Hermitian conjugate of r < p < q)
            // q < p < r (1/4)(1-Z_q)(Z_{r-1,p+1})(X_p X_r + Y_p Y_r)
            // p < q < r (1/4)(1-Z_q)(Z_{r-1,p+1})(X_p X_r + Y_p Y_r)
            
            // This amounts to applying a PQ term, followed by same PQ term after a CNOT from q to the parity bit.
            if (Length(idxFermions) == 2) {
                let termPR0 = GeneratorIndex((idxTermType, [1.0]), idxFermions);
                _ApplyJordanWignerPQTerm_(termPR0, angle, new Qubit[0], qubits);
            }
            else {
                
                if (idxFermions[0] < qubitQidx and qubitQidx < idxFermions[3]) {
                    let termPR1 = GeneratorIndex((idxTermType, [1.0]), [idxFermions[0], idxFermions[3] - 1]);
                    _ApplyJordanWignerPQTerm_(termPR1, angle, new Qubit[0], Excluding([qubitQidx], qubits));
                }
                else {
                    let termPR1 = GeneratorIndex((idxTermType, [1.0]), [0, idxFermions[3] - idxFermions[0]]);
                    _ApplyJordanWignerPQTerm_(termPR1, angle, [qubits[qubitQidx]], qubits[idxFermions[0] .. idxFermions[3]]);
                }
            }
        }
        
        adjoint invert;
        controlled distribute;
        controlled adjoint distribute;
    }
    
    
    /// # Summary
    /// Applies time-evolution by a PQRS term described by a `GeneratorIndex`.
    ///
    /// # Input
    /// ## term
    /// `GeneratorIndex` representing a PQRS term.
    /// ## stepSize
    /// Duration of time-evolution.
    /// ## extraParityQubits
    /// Optional parity qubits that flip the sign of time-evolution.
    /// ## qubits
    /// Qubits of Hamiltonian.
    operation _ApplyJordanWigner0123Term_ (term : GeneratorIndex, stepSize : Double, qubits : Qubit[]) : Unit {
        
        body (...) {
            let ((idxTermType, v0123), idxFermions) = term!;
            let angle = stepSize;
            let qubitsPQ = Subarray(idxFermions[0 .. 1], qubits);
            let qubitsRS = Subarray(idxFermions[2 .. 3], qubits);
            let qubitsPQJW = qubits[idxFermions[0] + 1 .. idxFermions[1] - 1];
            let qubitsRSJW = qubits[idxFermions[2] + 1 .. idxFermions[3] - 1];
            let ops = [[PauliX, PauliX, PauliX, PauliX], [PauliX, PauliX, PauliY, PauliY], [PauliX, PauliY, PauliX, PauliY], [PauliY, PauliX, PauliX, PauliY], [PauliY, PauliY, PauliY, PauliY], [PauliY, PauliY, PauliX, PauliX], [PauliY, PauliX, PauliY, PauliX], [PauliX, PauliY, PauliY, PauliX]];
            
            for (idxOp in IndexRange(ops)) {
                
                if (IsNotZero(v0123[idxOp % 4])) {
                    Exp(ops[idxOp] + ConstantArray(Length(qubitsPQJW) + Length(qubitsRSJW), PauliZ), angle * v0123[idxOp % 4], ((qubitsPQ + qubitsRS) + qubitsPQJW) + qubitsRSJW);
                }
            }
        }
        
        adjoint invert;
        controlled distribute;
        controlled adjoint distribute;
    }
    
    
    /// # Summary
    /// Converts a Hamiltonian described by `JWOptimizedHTerms`
    /// to a `GeneratorSystem` expressed in terms of the
    /// `GeneratorIndex` convention defined in this file.
    ///
    /// # Input
    /// ## data
    /// Description of Hamiltonian in `JWOptimizedHTerms` format.
    ///
    /// # Output
    /// Representation of Hamiltonian as `GeneratorSystem`.
    function JordanWignerGeneratorSystem (data : JWOptimizedHTerms) : GeneratorSystem {
        
        let (ZData, ZZData, PQandPQQRData, h0123Data) = data!;
        let ZGenSys = HTermsToGenSys(ZData, [0]);
        let ZZGenSys = HTermsToGenSys(ZZData, [1]);
        let PQandPQQRGenSys = HTermsToGenSys(PQandPQQRData, [2]);
        let h0123GenSys = HTermsToGenSys(h0123Data, [3]);
        return SumGeneratorSystems([ZGenSys, ZZGenSys, PQandPQQRGenSys, h0123GenSys]);
    }
    
    
    /// # Summary
    /// Represents a dynamical generator as a set of simulatable gates and an
    /// expansion in the JordanWigner basis.
    ///
    /// See [Dynamical Generator Modeling](../libraries/data-structures#dynamical-generator-modeling)
    /// for more details.
    ///
    /// # Input
    /// ## generatorIndex
    /// A generator index to be represented as unitary evolution in the JordanWigner.
    /// ## stepSize
    /// A multiplier on the duration of time-evolution by the term referenced
    /// in `generatorIndex`.
    /// ## qubits
    /// Register acted upon by time-evolution operator.
    operation JordanWignerFermionImpl (generatorIndex : GeneratorIndex, stepSize : Double, qubits : Qubit[]) : Unit {
        
        body (...) {
            let ((idxTermType, idxDoubles), idxFermions) = generatorIndex!;
            let termType = idxTermType[0];
            
            if (termType == 0) {
                _ApplyJordanWignerZTerm_(generatorIndex, stepSize, qubits);
            }
            elif (termType == 1) {
                _ApplyJordanWignerZZTerm_(generatorIndex, stepSize, qubits);
            }
            elif (termType == 2) {
                _ApplyJordanWignerPQandPQQRTerm_(generatorIndex, stepSize, qubits);
            }
            elif (termType == 3) {
                _ApplyJordanWigner0123Term_(generatorIndex, stepSize, qubits);
            }
        }
        
        adjoint invert;
        controlled distribute;
        controlled adjoint distribute;
    }
    
    
    /// # Summary
    /// Represents a dynamical generator as a set of simulatable gates and an
    /// expansion in the JordanWigner basis.
    ///
    /// # Input
    /// ## generatorIndex
    /// A generator index to be represented as unitary evolution in the JordanWigner.
    ///
    /// # Output
    /// An `EvolutionUnitary` representing time-evolution by the term
    /// referenced in `generatorIndex.
    function JordanWignerFermionFunction (generatorIndex : GeneratorIndex) : EvolutionUnitary {
        
        return EvolutionUnitary(JordanWignerFermionImpl(generatorIndex, _, _));
    }
    
    
    /// # Summary
    /// Represents a dynamical generator as a set of simulatable gates and an
    /// expansion in the JordanWigner basis.
    ///
    /// # Output
    /// An `EvolutionSet` that maps a `GeneratorIndex` for the JordanWigner basis to
    /// an `EvolutionUnitary.
    function JordanWignerFermionEvolutionSet () : EvolutionSet {
        
        return EvolutionSet(JordanWignerFermionFunction(_));
    }
    
}

