// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Chemistry.Tests {
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Diagnostics;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Preparation;
    open Microsoft.Quantum.Arithmetic;
    open Microsoft.Quantum.Math;
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Chemistry.JordanWigner;
    
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
        
        mutable product = [false, size = nFermions];
        mutable expected = [false, size = nFermions];
        for idx in 0..nFermions - 1 {
            set product w/= idx <- expectedBitString[idx] == bitString[idx] ? true | false;
            set expected w/= idx <- true;
        }
        
        AllEqualityFactB (product, expected, "Bit strings not equal");
    }

    function _JordanWignerClusterOperatorPQRSTermSignsTestHelper(idx : Int) : (Int[], Double[], Double){
        let cases = [
        ([1,2,3,4],[1.0,-1.0,-1.0,-1.0,1.0,1.0,1.0,-1.0],1.0),
        ([2,1,4,3],[1.0,-1.0,-1.0,-1.0,1.0,1.0,1.0,-1.0],1.0),
        ([3,4,1,2],[1.0,-1.0,-1.0,-1.0,1.0,1.0,1.0,-1.0],-1.0),
        ([2,1,3,4],[1.0,-1.0,-1.0,-1.0,1.0,1.0,1.0,-1.0],-1.0),
        ([1,3,2,4],[-1.0,-1.0,-1.0,1.0,-1.0,1.0,1.0,1.0],1.0),
        ([4,2,3,1],[-1.0,-1.0,-1.0,1.0,-1.0,1.0,1.0,1.0],-1.0),
        ([1,4,2,3],[1.0,1.0,-1.0,1.0,-1.0,1.0,-1.0,-1.0],1.0),
        ([2,3,4,1],[1.0,1.0,-1.0,1.0,-1.0,1.0,-1.0,-1.0],1.0)
        ];
        return cases[idx];
    }

    function _JordanWignerClusterOperatorPQRSTermSignsTest() : Unit{
        for idx in 0..7 {
            let (testCase, expectedSigns, expectedGlobalSign) = _JordanWignerClusterOperatorPQRSTermSignsTestHelper(idx);
            let (sortedIndices, signs, globalSign) = _JordanWignerClusterOperatorPQRSTermSigns(testCase);

            let p = sortedIndices[0];
            let q = sortedIndices[1];
            let r = sortedIndices[2];
            let s = sortedIndices[3];

            Fact(p < q and q < r and r < s, "Expected p < q < r < s.");
            NearEqualityFactD(globalSign, expectedGlobalSign);
            for (actual, expected) in Zipped(signs, expectedSigns)  {
                NearEqualityFactD(actual, expected);
            }
        }
    }

    function _DoublesToComplexPolar(input: Double[]) : ComplexPolar[]{
        mutable arr = [ComplexPolar(0.0, 0.0), size = Length(input)];
        for idx in 0..Length(input)-1 {
            set arr w/= idx <- ComplexAsComplexPolar(Complex(input[idx], 0.));
        }
        return arr;
    }

    operation _JordanWignerUCCTermTestHelper(nQubits: Int, excitations: Int[], term: JordanWignerInputState[], result: Double[]) : Unit{
        use qubits = Qubit[nQubits];
        for idx in excitations {
            X(qubits[idx]);
        }
        PrepareUnitaryCoupledClusterState(NoOp, term, 1.0, qubits);
        DumpRegister((), qubits);
        Adjoint PrepareArbitraryStateCP(_DoublesToComplexPolar(result), LittleEndian(qubits));
        AssertAllZeroWithinTolerance(qubits, 1e-5);
        ResetAll(qubits);
    }

    operation JordanWignerUCCSTermTest() : Unit{
        // test using Exp(2.0 (a^\dag_1 a_3 - h.c.)) 
        let term0 = [JordanWignerInputState((2.0,0.0), [1,3])];
        let state0 = [0.,0.,-0.416147,0.,0.,0.,0.,0.,-0.909297,0.,0.,0.,0.,0.,0.,0.];
        _JordanWignerUCCTermTestHelper(4, [1], term0, state0);

        // test using Exp(2.0 (a^\dag_3 a_1 - h.c.)) 
        let term1 = [JordanWignerInputState((2.0,0.0), [3,1])];
        let state1 = [0.,0.,-0.416147,0.,0.,0.,0.,0.,0.909297,0.,0.,0.,0.,0.,0.,0.];
        _JordanWignerUCCTermTestHelper(4, [1], term1, state1);
    }

    operation JordanWignerUCCDTermPQRSTest() : Unit{
        // test using Exp(2.0 (a^\dag_0 a^\dag_1 a_3 a_4 - h.c.)) 
        let term0 = [JordanWignerInputState((2.0,0.0), [0,1,2,4])];
        let state0 = [0.,0.,0.,-0.416147,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.909297,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.];
        _JordanWignerUCCTermTestHelper(5, [0,1], term0, state0);

        // test using Exp(2.0 (a^\dag_0 a^\dag_1 a_3 a_4 - h.c.)) 
        let term1 = [JordanWignerInputState((2.0,0.0), [0,1,2,4])];
        let state1 = [0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,-0.416147,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,-0.909297,0.,0.,0.];
        _JordanWignerUCCTermTestHelper(5, [0,1,3], term1, state1);

        // test using Exp(2.0 (a^\dag_1 a^\dag_0 a_2 a_4 - h.c.)) 
        let term2 = [JordanWignerInputState((2.0,0.0), [1,0,2,4])];
        let state2 = [0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,-0.416147,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.909297,0.,0.,0.];
        _JordanWignerUCCTermTestHelper(5, [0,1,3], term2, state2);

        // test using Exp(2.0 (a^\dag_1 a^\dag_0 a_2 a_4 - h.c.)) 
        let term3 = [JordanWignerInputState((-2.0,0.0), [4,2,1,0])];
        let state3 = [0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,-0.416147,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.909297,0.,0.,0.];
        _JordanWignerUCCTermTestHelper(5, [0,1,3], term2, state2);
    }

    operation JordanWignerUCCDTermPRQSTest() : Unit {
        let term0 = [JordanWignerInputState((2.0,0.0), [2,0,4,1])];
        let state0 = [0.,0.,0.,0.,0.,-0.416147,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.909297,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.];
        _JordanWignerUCCTermTestHelper(5, [0,2], term0, state0);

        let term1 = [JordanWignerInputState((2.0,0.0), [2,0,4,1])];
        let state1 = [0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,-0.416147,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,-0.909297,0.,0.,0.,0.,0.];
        _JordanWignerUCCTermTestHelper(5, [0,2,3], term1, state1);
    }

    operation JordanWignerUCCDTermPRSQTest() : Unit {
        let term3 = [JordanWignerInputState((2.0,0.0), [0,4,2,3])];
        let state3 = [0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.909297,0.,0.,0.,0.,-0.416147,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.];
        _JordanWignerUCCTermTestHelper(5, [0,4], term3, state3);
    }

}
