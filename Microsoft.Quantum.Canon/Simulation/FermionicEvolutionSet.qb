// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

// FIXME: This file is still a work in progress, and should not yet be included
//        in builds.
namespace Microsoft.Quantum.Canon {
    open Microsoft.Quantum.Primitive

    // Convention for GeneratorIndex = ((Int[],Double[]), Int[])
    // We index annihilation operators as 0 and creation operators as 1.
    // For instance, the Hermitian Fermionic string a_p^\dag a_q^\dag a_r a_s + Hermitian Conjugate
    // would be indexed as a = [1,1,0,0].
    // We do not use the Double[] for indexing Pauli strings e.g. b = new Double[0]
    // We index the Fermionics that Fermionic strings act on with arrays of integers e.g. f = [p,q,e,s] 
    // An example of a Fermionic string GeneratorIndex is thus ((a,b), q)

    // Consider the Hamiltonian H = 0.1 (a_2^\dag) + 0.2 (a_0^\dag a_2) + h.c.
    // Its GeneratorTerms are (([1],[]),[2]), 0.1),  (([1,0],[]),[0,2]), 0.2).

    /// TODO Implement liquid PP | PQ | PQQP | PQQR | PQRS

    /// summary:
    ///     Represents a dynamical generator as a set of simulatable gates and
    ///     an expansion in terms of that basis.
    ///     Last parameter for number of terms
    ///     Converts creation and annihilation operators to Pauli string

    function JordanWignerZString(idxFermionMin : Int, idxFermionMax : Int, coefficient: Double) : (Int[], Int[])
    {
        let nQubits = idxFermionMax - idxFermionMin + 1
        mutable idxPauliString = new Int[nQubits]
        mutable idxQubits = new Int[nQubits]
        for(idxQubit in idxQubitMin..idxQubitMax){
            set idxPauliString[idxQubit] = 3
            set idxQubits[idxQubit] = idxQubit
        }
        return (idxPauliString, idxQubits)
    }

    function JordanWignerZString(idxFermionMin : Int, idxFermionMax : Int, coefficient: Double) : (Int[], Int[])
    {
        let nQubits = idxFermionMax - idxFermionMin + 1
        mutable idxPauliString = new Int[nQubits]
        mutable idxQubits = new Int[nQubits]
        for(idxQubit in idxQubitMin..idxQubitMax){
            set idxPauliString[idxQubit] = 3
            set idxQubits[idxQubit] = idxQubit
        }
        return (idxPauliString, idxQubits)
    }

    /// # References
    /// - https://arxiv.org/pdf/1001.3855.pdf
    function JordanWigner(generatorIndex : GeneratorIndex) : GeneratorIndex[] {
        let (idxGen, idxFermions) = generatorIndex
        let (idxCreationAnnihilation, unused) = idxGen
        let nOperators = Length(idxCreationAnnihilation)


        if(nOperators == 2){
            /// Number operator p^d p
            if( idxFermions[0] == idxFermions[1] &&             
                idxCreationAnnihilation[0] == 1 && 
                idxCreationAnnihilation[1] == 0 &&)
            {
                let idxQubit = idxFermions[0]
                R1(theta, qubits[idxQubit])
            }
            elif(idxFermions[0] != idxFermions[1] &&             
                idxCreationAnnihilation[0] == 1 && 
                idxCreationAnnihilation[1] == 0 &&)
            {
                let idxQubitP = Min(idxFermions)
                let idxQubitQ = Max(idxFermions)
                let qubitP = qubits[idxQubitP]
                let qubitQ = qubits[idxQubitQ]
                H(qubitP)
                H(qubitQ)
                R1(theta, qubits[idxQubit])
            }
        }

    }

    function JordanWigner(generatorIndex : GeneratorIndex) : GeneratorIndex[] {
        let (idxGen, idxFermions) = generatorIndex
        let (idxCreationAnnihilation, unused) = idxGen
        let nOperators = Length(idxCreationAnnihilation)
        let e = new Double[0]
        // Single-site a_k^\dag + a_k
        if(nOperators == 1){
            mutable idxPauliString = new Int[idxFermions[0] + 1]
            mutable idxQubits = new Int[idxFermions[0] + 1]
            for(idxQubit in 0..idxFermions[0] - 1){
                set idxPauliString[idxQubit] = 3
                set idxQubits[idxQubit] = idxQubit
            }
            set idxPauliString[idxFermions[0]] = 1
            set idxQubits[idxFermions[0]] = idxFermions[0]

            let generatorIndex0 = GeneratorIndex((idxPauliString, e), idxQubits)

            return [generatorIndex0]
        }
        // PP | PQ a_p^\dag a_q + a_q^\dag a_p
        elif(nOperators == 2){
            let idxQubitMax = Max(idxFermions)
            let idxQubitMin = Min(idxFermions)
            let nQubits = idxQubitMax - idxQubitMin + 1
            mutable idxPauliString = new Int[nQubits]
            mutable idxQubits = new Int[nQubits]
            for(idxQubit in idxQubitMin + 1..idxQubitMax - 1){
                set idxPauliString[idxQubit] = 3
                set idxQubits[idxQubit] = idxQubit
            }
            set idxPauliString[idxQubitMin] = 1
            set idxPauliString[idxQubitMax] = 1
            set idxQubits[idxQubitMin] = idxQubitMin
            set idxQubits[idxQubitMax] = idxQubitMax

            let generatorIndex0 = GeneratorIndex((idxPauliString, e), idxQubits)

            set idxPauliString[idxQubitMin] = 2
            set idxPauliString[idxQubitMax] = 2

            let generatorIndex1 = GeneratorIndex((idxPauliString, e), idxQubits)

            return [generatorIndex0;generatorIndex1]
        }
        //PQQP | PQQR | PQRS
        elif(nOperators == 4){
            // p^\dag p q^\dag q
            if( idxFermions[0] == idxFermions[1] && 
                idxFermions[2] == idxFermions[3] && 
                idxCreationAnnihilation[0] == 1 && 
                idxCreationAnnihilation[1] == 0 && 
                idxCreationAnnihilation[2] == 1 &&
                idxCreationAnnihilation[3] == 0 ){
                
            }
            // p^d q^d q p
            elif( idxFermions[0] == idxFermions[3] && 
                idxFermions[1] == idxFermions[2] && 
                idxCreationAnnihilation[0] == 1 && 
                idxCreationAnnihilation[1] == 1 && 
                idxCreationAnnihilation[2] == 9 &&
                idxCreationAnnihilation[3] == 0 ){
                
                let idxP = idxFermions[0]
                let idxQ = idxFermions[1]

                [GeneratorIndex(([0], e), [idxP]);
                GeneratorIndex(([0], e), [idxP]);]

                let idxQubitMax = Max(idxFermions)
                let idxQubitMin = Min(idxFermions)
                let nQubits = idxQubitMax - idxQubitMin + 1
                mutable idxPauliString = new Int[nQubits]
                mutable idxQubits = new Int[nQubits]
                for(idxQubit in idxQubitMin + 1..idxQubitMax - 1){
                    set idxPauliString[idxQubit] = 3
                    set idxQubits[idxQubit] = idxQubit
                }

            }
            

        }
        else{
            // TODO for 3 terms and 4 terms
            let generatorIndex0 = new GeneratorIndex[0]
            return generatorIndex0
        }
    }

    function SelectFermionicEncoding(fermionicEncoding : Int, generatorIndex : GeneratorIndex) : GeneratorIndex[]
    {
        //C hoose encoding from Fermionic operators to qubits. Defaults to Jordan-Wigner
        if(fermionicEncoding == 0){
            return JordanWigner(generatorIndex)
        }
        else{
            return JordanWigner(generatorIndex)
        }
    }

    // This calls EGPauli.qb to implement the Jordan-Wigner strings.
    operation FermionicEvolutionSetImpl(fermionicEncoding : Int, generatorIndex : GeneratorIndex, delta : Double, qubits: Qubit[]) : ()
    {
        body{
            let generatorIndexArr = SelectFermionicEncoding(fermionicEncoding, generatorIndex)
            let nGeneratorIndex = Length(generatorIndexArr)
            for(idxGeneratorIndex in 0..nGeneratorIndex-1){
                (PauliEvolutionSet(generatorIndexArr[idxGeneratorIndex]))(delta,qubits)
            }
        }
        adjoint auto
        controlled auto
        controlled adjoint auto
    }

    function FermionicEvolutionSetImpl2(fermionicEncoding : Int, generatorIndex : GeneratorIndex) : EvolutionUnitary
    {
        return EvolutionUnitary(FermionicEvolutionSetImpl(fermionicEncoding, generatorIndex, _, _))
    }

    function FermionicEvolutionSet(fermionicEncoding : Int) : (GeneratorIndex -> EvolutionUnitary)
    {
        return FermionicEvolutionSetImpl2(fermionicEncoding,_)
    }

}
