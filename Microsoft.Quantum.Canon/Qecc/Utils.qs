// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {
    open Microsoft.Quantum.Primitive;

    /// # Summary
    /// Given an array of results, represents the array by a single
    /// integer, with the 0th (leftmost) entry in the array being mapped
    /// the least significant bit. Thus, [One; Zero] is represented by
    /// 1 and [Zero; One] by 2.
    function ResultAsInt(results : Result[])  : Int
    {
        mutable n = 0;

        for (idxResult in 0..(Length(results) - 1)) {
            if (results[idxResult] == One) {
                set n = n + 2 ^ idxResult;
            }
        }

        return n;
    }

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
	/// This does not checks if the given set of generators are commuting.
	/// If they are not commuting, the result of measurements may be arbitrary.
    operation  MeasureStabilizerGenerators(stabilizerGroup : Pauli[][],  logicalRegister : LogicalRegister, gadget : ((Pauli[], Qubit[]) => Result))  : Syndrome
    {
        body {
            let results = MeasurePaulis(stabilizerGroup, logicalRegister, gadget);
            return Syndrome(results);
        }
    }

    operation Recover( code : QECC,  fn : RecoveryFn,  logicalRegister : LogicalRegister)  : ()
    {
        body {
            let (encode, decode, syndMeas) = code;
            let syndrome = syndMeas(logicalRegister);
            let recoveryOp = fn(syndrome);
            ApplyPauli(recoveryOp, logicalRegister);
        }
    }

    operation RecoverCSS( code : CSS,  fnX : RecoveryFn,  fnZ : RecoveryFn,  logicalRegister : LogicalRegister)  : ()
    {
        body {
            let (encode, decode, syndMeasX, syndMeasZ) = code;
            let syndromeX = syndMeasX(logicalRegister);
            let recoveryOpX = fnX(syndromeX);
            Message($"X: {syndromeX} → {recoveryOpX}");
            ApplyPauli(recoveryOpX, logicalRegister);
            let syndromeZ = syndMeasZ(logicalRegister);
            let recoveryOpZ = fnZ(syndromeZ);
            Message($"Z: {syndromeZ} → {recoveryOpZ}");
            ApplyPauli(recoveryOpZ, logicalRegister);
        }
    }

    function TableLookupRecoveryImpl(table : Pauli[][],  syndrome : Syndrome)  : Pauli[]
    {
        return table[ResultAsInt(syndrome)];
    }

    function  TableLookupRecovery(table : Pauli[][])  : RecoveryFn
    {
        return RecoveryFn(TableLookupRecoveryImpl(table, _));
    }

}
