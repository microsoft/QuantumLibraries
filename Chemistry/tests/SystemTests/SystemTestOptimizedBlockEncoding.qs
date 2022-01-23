// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace SystemTestsOptimizedBlockEncoding {

    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Convert;
    open Microsoft.Quantum.Chemistry.JordanWigner;

    //////////////////////////////////////////////////////////////////////////
    // Using OptimizedBlockEncoding() ////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////

    /// # Summary
    /// We can now use Canon's phase estimation algorithms to
    /// learn the ground state energy using the above simulation.
    operation RunOptimizedBlockEncoding (nSpinOrbitals : Int, data : JWOptimizedHTerms, targetError : Double) : Unit {
        let ((nCtrlRegisterQubits, nTargetRegisterQubits), (oneNorm, blockEncodingReflection)) = _JordanWignerOptimizedBlockEncoding_(targetError, data, nSpinOrbitals);
        let op = blockEncodingReflection!!;

        use ctrlQubits = Qubit[nCtrlRegisterQubits];
        use targetQubits = Qubit[nTargetRegisterQubits];
        op(ctrlQubits, targetQubits);
    }

}
