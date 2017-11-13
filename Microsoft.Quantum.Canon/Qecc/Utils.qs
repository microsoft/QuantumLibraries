// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {

    /// <summary>
    ///     Given an array of results, represents the array by a single
    ///     integer, with the 0th (leftmost) entry in the array being mapped
    ///     the least significant bit. Thus, [One; Zero] is represented by
    ///     1 and [Zero; One] by 2.
    /// </summary>
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
            ApplyPauli(recoveryOpX, logicalRegister);
            let syndromeZ = syndMeasZ(logicalRegister);
            let recoveryOpZ = fnZ(syndromeZ);
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
