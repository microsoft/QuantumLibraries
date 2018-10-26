// Copyright (c) Microsoft Corporation. All rights reserved. Licensed under the
// Microsoft Software License Terms for Microsoft Quantum Simulation Library (Preview).
// See LICENSE.md in the project root for license information.
namespace SystemTestsOptimizedBlockEncoding {
    
    open Microsoft.Quantum.Primitive;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Extensions.Math;
    open Microsoft.Quantum.Extensions.Convert;
    open Microsoft.Quantum.Chemistry.JordanWigner;  
    
    //////////////////////////////////////////////////////////////////////////
    // Using OptimizedBlockEncoding() ////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////
    
    /// # Summary
    /// We can now use Canon's phase estimation algorithms to
    /// learn the ground state energy using the above simulation.
    operation RunOptimizedBlockEncoding (nSpinOrbitals : Int, data : JWOptimizedHTerms, targetError : Double) : Unit {
        
        let ((nCtrlRegisterQubits, nTargetRegisterQubits), (oneNorm, blockEncodingReflection)) = _JordanWignerOptimizedBlockEncoding_(targetError, data, nSpinOrbitals);
        let nQubits = nCtrlRegisterQubits + nTargetRegisterQubits;
        let op = blockEncodingReflection!!;
        using (ctrlQubits = Qubit[nCtrlRegisterQubits]) {
            
            using (targetQubits = Qubit[nTargetRegisterQubits]) {
                op(ctrlQubits, targetQubits);
                //ResetAll(ctrlQubits + targetQubits);
            }
        }
    }
    
}


