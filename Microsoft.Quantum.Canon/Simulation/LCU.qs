// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

// FIXME: This file is still a work in progress, and should not yet be included
//        in builds.
namespace Microsoft.Quantum.Canon {
    // NB: Needed for H and X in Qubitization.
    open Microsoft.Quantum.Primitive
    open Microsoft.Quantum.Extensions.Math;

    //Intent is to provide generic construction of LCU & qubitization technique
    //Example implementation in terms of PauliStrings to be given, which will also cover Fermionic encodings
    //Qubitization approach
    newtype UnitaryOperator = (Qubit[] => () : Adjoint, Controlled)
    newtype OracleLCUSelector = ((Qubit[], Qubit[]) => () : Adjoint, Controlled)
    newtype OracleLCUStatePrep = (Qubit[] => () : Adjoint, Controlled)

    //Need to write op for arbitrary state preparation of |0> -> a_1|0> + a_2|1> + a_3|2> +...
    operation LCUStatePrepImpl ( amplitude: Double[], qubits : Qubit[]) : (){
        body{



        }

        adjoint auto
        controlled auto
    }
    function LCUStatePrep(amplitude : Double[]) : OracleLCUStatePrep
    {
        return LCUStatePrepImpl(amplitude, _)
    }



    //Selector is multiply-controlled unitary multiplexor
    operation LCUSelectorImpl ( unitaryOperators :  UnitaryOperator[] , qubitsAncilla : Qubit[], qubitsSystem : Qubit[]  ): (){
        body{
            let nUnitaries = Length(unitaryOperators)
            for(idxUnitary in 0..nUnitaries){
                ControlledOnInt(idxUnitary , unitaryOperators[idxUnitary])(qubitsAncilla, qubitsSystem)
            }
        }
    }
    function LCUSelector ( unitaryOperators :  UnitaryOperator[] ) : OracleLCUSelector {
        return LCUSelectorImpl(unitaryOperators, _, _)
    }

    //Example implementation of unitaryOperator given description as PauliString
    operation UnitaryOperatorFromPauliStringImpl ( pauliString: PauliString, qubits : Qubit[]) : (){
        body{
            let (paulis, qubitIndices) = pauliString
            ApplyPauli(paulis, TEMPQubitSlice(qubitIndices, qubits))
        }
    }
    function UnitaryOperatorFromPauliString( pauliString: PauliString) : UnitaryOperator
    {
        return UnitaryOperatorFromPauliStringImpl ( pauliString,_)
    }

    //Selector for Paulistrings
    operation LCUSelectorPauliImpl ( pauliStrings : PauliString[] , qubitsAncilla : Qubit[], qubitsSystem : Qubit[] ) : ()
    {
        body{
            let nUnitaries = Length(pauliString)
            mutable unitaries = new UnitaryOperator[nUnitaries]
            for(idx in 0.. nUnitaries-1){
                unitaries[idx] = UnitaryOperatorFromPauliString(pauliStrings[idx])
            }
            LCUSelector(unitaries , qubitsAncilla, qubitsSystem)
        }
    }
    function LCUSelectorPauli ( pauliStrings : PauliString[] ) : OracleLCUSelector
    {
        return LCUSelectorPauliImpl(pauliStrings, _, _)
    }

    //Qubitization creates quantum walk -- eigenphases related to eigenvalues of Hamiltonian like ArcSin[Hamiltonian/Normalization]
    //Note that 1 qubit less can be used if PauliStrings rather than arbitrary unitaries areused
    operation Qubitization( oracleLCUState : OracleLCUStatePrep, oracleLCUSelect: OracleLCUSelector , qubitQubitization : Qubit, qubitsAncilla : Qubit[], qubitsSystem : Qubit[] ) : ()
    {
        body {
            H(qubitQubitization)
            oracleLCUState(qubitsAncilla)
            (Controlled oracleLCUSelect)([qubitQubitization], qubitsAncilla, qubitsSystem)
            X(qubitQubitization)
            (Controlled Adjoint oracleLCUSelect)([qubitQubitization], qubitsAncilla, qubitsSystem)
            (Adjoint oracleLCUState)(qubitsAncilla)
            H(qubitQubitization)
            ReflectionAllZero(PI(), [qubitQubitization] + qubitsAncilla)
        }
    }



    //Example find g.s. energy of H2
    operation ExampleH2ByLCU(qubitSystem : Qubit[]) : (){
        body {
            let idxBond = 1
            let (coeffs, paulis) = HamH2(idxBond)

            //apply square root to all coeffs
            mutable coeffsSqrt = coeffs 
            nCoeffs = Length(coeffs)
            let oracleLCUState = LCUStatePrep (coeffsSqrt)
            let oracleLCUSelector = LCUSelectorPauli ( paulis)

            using(qubitsExtra = Qubit[1+ nCoeffs]){

                let qubiterate =  Qubitization(oracleLCUState, oracleLCUSelector, qubitsExtra[0], qubitsExtra[1..nCoeffs - 1], qubitsSystem)
                //Perform QPE on qubiterate
            }
        }
    }

}
