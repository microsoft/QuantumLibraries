// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.ErrorCorrection {
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Measurement;
    open Microsoft.Quantum.Convert;

    /// # Summary
    /// Measures the given set of generators of a stabilizer group.
    /// # Input
    /// ## stabilizerGroup
    /// An array of multiqubit Pauli operators.
    /// For example, `stabilizerGroup[0]` is a list of single-qubit Pauli matrices,
    /// the tensor product of which is a stabilizer generator.
    /// ## logicalRegister
    /// An array of qubits where the stabilizer code is defined.
    /// ## gadget
    /// An operation that specifies how to measure a multiqubit Pauli operator.
    /// # Output
    /// The result of measurements.
    /// # Remarks
    /// This does not check if the given set of generators are commuting.
    /// If they are not commuting, the result of measurements may be arbitrary.
    operation MeasureStabilizerGenerators (stabilizerGroup : Pauli[][], logicalRegister : LogicalRegister, gadget : ((Pauli[], Qubit[]) => Result)) : Syndrome
    {
        let results = MeasurePaulis(stabilizerGroup, logicalRegister!, gadget);
        return Syndrome(results);
    }
    
    
    /// # Summary
    /// Performs a single round of error correction by a quantum code
    /// described by a `QECC` type.
    ///
    /// # Input
    /// ## code
    /// A quantum error-correcting code packaged as a `QECC` type describes
    /// the encoding and decoding of quantum data, and how error syndromes
    /// are to be measured.
    /// ## fn
    /// A `RecoveryFn` that maps each error syndrome to the `Pauli[]` operations
    /// that correct the detected error.
    /// ## logicalRegister
    /// An array of qubits where the stabilizer code is defined.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.QECC
    /// - Microsoft.Quantum.Canon.RecoveryFn
    /// - Microsoft.Quantum.Canon.LogicalRegister
    operation Recover (code : QECC, fn : RecoveryFn, logicalRegister : LogicalRegister) : Unit
    {
        let (encode, decode, syndMeas) = code!;
        let syndrome = syndMeas!(logicalRegister);
        let recoveryOp = fn!(syndrome);
        ApplyPauli(recoveryOp, logicalRegister!);
    }
    
    
    /// # Summary
    /// Performs a single round of error correction by a quantum code
    /// described by a `CSS` type.
    ///
    /// # Input
    /// ## code
    /// A quantum CSS error-correcting code packaged as a `CSS` type describes
    /// the encoding and decoding of quantum data, and how error syndromes
    /// are to be measured.
    /// ## fnX
    /// A `RecoveryFn` that maps each $X$ error syndrome to the `Pauli[]` operations
    /// that correct the detected error.
    /// ## fnZ
    /// A `RecoveryFn` that maps each $Z$ error syndrome to the `Pauli[]` operations
    /// that correct the detected error.
    /// ## logicalRegister
    /// An array of qubits where the stabilizer code is defined.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.CSS
    /// - Microsoft.Quantum.Canon.RecoveryFn
    /// - Microsoft.Quantum.Canon.LogicalRegister
    operation RecoverCSS (code : CSS, fnX : RecoveryFn, fnZ : RecoveryFn, logicalRegister : LogicalRegister) : Unit
    {
        let (encode, decode, syndMeasX, syndMeasZ) = code!;
        let syndromeX = syndMeasX!(logicalRegister);
        let recoveryOpX = fnX!(syndromeX);
        Message($"X: {syndromeX} → {recoveryOpX}");
        ApplyPauli(recoveryOpX, logicalRegister!);
        let syndromeZ = syndMeasZ!(logicalRegister);
        let recoveryOpZ = fnZ!(syndromeZ);
        Message($"Z: {syndromeZ} → {recoveryOpZ}");
        ApplyPauli(recoveryOpZ, logicalRegister!);
    }
    
    
    function TableLookupRecoveryImpl (table : Pauli[][], syndrome : Syndrome) : Pauli[]
    {
        return table[ResultArrayAsInt(syndrome!)];
    }
    
    
    /// # Summary
    /// For a given table of Pauli operations on a given register of qubits, this function
    /// returns an object of type `RecoveryFn` which contains all information needed to
    /// perform a table-lookup decoding with respect to the given array of Pauli operations.
    ///
    /// # Input
    /// ## table
    /// Table of Pauli operations that operate on a given qubit register
    ///
    /// # Output
    /// An object of type `RecoveryFn`, i.e., a map `Syndrome -> Pauli[]` that associates
    /// with a given syndrome array a corresponding Pauli correction operation.
    function TableLookupRecovery (table : Pauli[][]) : RecoveryFn
    {
        return RecoveryFn(TableLookupRecoveryImpl(table, _));
    }
    
}


