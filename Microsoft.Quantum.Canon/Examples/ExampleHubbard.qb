// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Examples.Hubbard {
    open Microsoft.Quantum.Canon

    // FIXME: TO be rewritten to match new library

    // Hubbard model example using the Hamiltonian representation library
    // Probably need to include spin index?
    // H = (-t \Sum_j  a_j^\dag a_{j+1} + h.c. ) + u \Sum_j  a_j^\dag a_j a_j^\dag a_j
    // Ising model using the Hamiltonian representation library
    // idxHamiltonian is in [0, 3 * nsites - 2]

    function HubbardGeneratorSystemImpl(nsites : Int, t : Double, u: Double, idxHamiltonian : Int) : GeneratorTerm
    {
        // when idxHamiltonian is in [0, nsites - 1], apply repulsion term "u"
        // when idxHamiltonian is in [nsites, 2 * nsites - 2], apply hopping term "t"
            let e = new Double[0]

            if(idxHamiltonian <= nsites - 1){
                let idxFermionicString = [1,0,1,0]
                let idxFermions = [idxHamiltonian, idxHamiltonian, idxHamiltonian, idxHamiltonian]
                let coeff = u
            }
            if(idxHamiltonian <= 2 * nsites - 2){
                let idxFermionicString = [1,0]
                let idxFermions = [idxHamiltonian % nsites, (idxHamiltonian % nsites) + 1]
                let coeff = -1.0 * t
            }
            else{
                let idxFermionicString = new Int[0]
                let idxFermions = new Int[0]
                let coeff = 0.0
            }
            let generatorIndex = GeneratorIndex((idxFermionicString,e),idxFermions)
            let generatorTerm = GeneratorTerm(generatorIndex, coeff)
            return generatorTerm
    }

    function HubbardGeneratorSystem(nsites : Int, t : Double, u: Double) : GeneratorSystem 
    {
        return GeneratorSystem(HubbardGeneratorSystemImpl(nsites, t, u, _))
    }

    functionHubbardEvolutionGenerator(nsites : Int, t : Double, u: Double) : EvolutionGenerator
    {
        let nTerms = 2 * nsites - 1
        let generatorSystem = HubbarGeneratorSystem(nsites, t, u)
        let fermionicEncoding = 0
        let evolutionSet = FermionicEvolutionSet(fermionicEncoding)
        return EvolutionGenerator(evolutionSet, generatorSystem, nTerms)
    }


    operation HubbardTrotterStep(nsites : Int, t : Double, u: Double, trotterOrder: Int, trotterStepSize: Double): (Qubit[] => ())
    {
        body{
            let evolutionGenerator = HubbardEvolutionGenerator(nsites, t, u)
            let (evolutionSet, generatorSystem, nTerms) = evolutionGenerator
            let op = ApplyEvolution(evolutionGenerator)
            (TrotterCA(op,nTerms,trotterOrder))(trotterStepSize, _)
    }

}
