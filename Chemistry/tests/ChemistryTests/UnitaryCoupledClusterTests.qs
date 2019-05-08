// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Chemistry.Tests {
    
    open Microsoft.Quantum.Primitive;
    open Microsoft.Quantum.Extensions.Testing;
    open Microsoft.Quantum.Extensions.Math;
    open Microsoft.Quantum.Extensions.Convert;
    open Microsoft.Quantum.Chemistry.JordanWigner; 
    open Microsoft.Quantum.Diagnostics;
    
    // Check that correct Pauli Z string is computed
    function _ComputeJordanWignerBitString_0Test() : Unit {
        let nFermions = 5;
        let fermionIndices = [0, 3];
        let expectedBitString = [false, true, true, false,false];
        let bitString = _ComputeJordanWignerBitString(nFermions, fermionIndices);
        AllEqualityFactB (bitString, expectedBitString, "Bit strings not equal");
    }

    function _ComputeJordanWignerBitString_1Test() : Unit {
        let nFermions = 7;
        let fermionIndices = [0, 4, 2, 6];
        let expectedBitString = [false, true, false, false, false, true, false];
        let bitString = _ComputeJordanWignerBitString(nFermions, fermionIndices);
        AllEqualityFactB (bitString, expectedBitString, "Bit strings not equal");
    }

    function _ComputeJordanWignerPauliZString_0Test() : Unit {
        let nFermions = 7;
        let fermionIndices = [0, 4, 2, 6];
        let expectedBitString = [PauliI, PauliZ, PauliI, PauliI, PauliI, PauliZ, PauliI];
        let bitString = _ComputeJordanWignerPauliZString(nFermions, fermionIndices);
        
        mutable product = new Bool[nFermions];
        mutable expected = new Bool[nFermions];
        for(idx in 0..nFermions - 1){
            set product w/= idx <- expectedBitString[idx] == bitString[idx] ? true | false;
            set expected w/= idx <- true;
        }
        
        AllEqualityFactB (product, expected, "Bit strings not equal");
    }
}


