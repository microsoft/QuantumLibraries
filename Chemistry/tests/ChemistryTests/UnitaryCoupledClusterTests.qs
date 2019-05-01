// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Chemistry.Tests {
    
    open Microsoft.Quantum.Primitive;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Extensions.Testing;
    open Microsoft.Quantum.Extensions.Math;
    open Microsoft.Quantum.Extensions.Convert;
    open Microsoft.Quantum.Chemistry.JordanWigner; 
    
    // Check that correct Pauli Z string is computed
    function _ComputeJordanWignerBitString_0Test() : Unit {
        let nFermions = 5;
        let fermionIndices = [0, 3];
        let expectedBitString = [false, true, true, false,false];
        let bitString = _ComputeJordanWignerBitString(nFermions, fermionIndices);
        AssertBoolArrayEqual (bitString, expectedBitString, "Bit strings not equal");
    }

    function _ComputeJordanWignerBitString_1Test() : Unit {
        let nFermions = 7;
        let fermionIndices = [0, 4, 2, 6];
        let expectedBitString = [false, true, false, false, false, true, false];
        let bitString = _ComputeJordanWignerBitString(nFermions, fermionIndices);
        AssertBoolArrayEqual (bitString, expectedBitString, "Bit strings not equal");
    }

    function _ComputeJordanWignerPauliZString_0Test() : Unit {
        let nFermions = 7;
        let fermionIndices = [0, 4, 2, 6];
        let expectedBitString = [PauliI, PauliZ, PauliI, PauliI, PauliI, PauliZ, PauliI];
        let bitString = _ComputeJordanWignerPauliZString(nFermions, fermionIndices);
        
        mutable product = new Bool[nFermions];
        mutable expected = new Bool[nFermions];
        for(idx in 0..nFermions - 1){
            set product[idx] = expectedBitString[idx] == bitString[idx] ? true | false;
            set expected[idx] = true;
        }
        
        AssertBoolArrayEqual (product, expected, "Bit strings not equal");
    }

    function _JordanWignerClusterOperatorPQRSTermSignsTestHelper(idx : Int) : (Int[], Double[], Double){
        mutable cases = new (Int[], Double[], Double)[8];
        set cases [0] = ([1,2,3,4],[1.0,-1.0,-1.0,-1.0,1.0,1.0,1.0,-1.0],1.0);
        set cases [1] = ([2,1,4,3],[1.0,-1.0,-1.0,-1.0,1.0,1.0,1.0,-1.0],1.0);
        set cases [2] = ([3,4,1,2],[1.0,-1.0,-1.0,-1.0,1.0,1.0,1.0,-1.0],-1.0);
        set cases [3] = ([2,1,3,4],[1.0,-1.0,-1.0,-1.0,1.0,1.0,1.0,-1.0],-1.0);
        set cases [4] = ([1,3,2,4],[-1.0,-1.0,-1.0,1.0,-1.0,1.0,1.0,1.0],1.0);
        set cases [5] = ([4,2,3,1],[-1.0,-1.0,-1.0,1.0,-1.0,1.0,1.0,1.0],-1.0);
        set cases [6] = ([1,4,2,3],[1.0,1.0,-1.0,1.0,-1.0,1.0,-1.0,-1.0],1.0);
        set cases [7] = ([2,3,4,1],[1.0,1.0,-1.0,1.0,-1.0,1.0,-1.0,-1.0],1.0);
        return cases[idx];
    }

    function _JordanWignerClusterOperatorPQRSTermSignsTest() : Unit{
        for(idx in 0..7){
            let (testCase, expectedDoubles, expectedGlobalSign) = _JordanWignerClusterOperatorPQRSTermSignsTestHelper(idx);
            let (sortedIndices, signs, globalsign) = _JordanWignerClusterOperatorPQRSTermSigns(testCase);

            let p = sortedIndices[0];
            let q = sortedIndices[1];
            let r = sortedIndices[2];
            let s = sortedIndices[3];

            AssertBoolEqual(true, p<q and q<r and r<s, "Expected p<q<r<s");
            AssertAlmostEqual(globalsign, expectedGlobalSign);
            for(signIdx in 0..Length(signs)-1){
                AssertAlmostEqual(signs[signIdx], expectedDoubles[signIdx]);
            }
        }
    }
}



