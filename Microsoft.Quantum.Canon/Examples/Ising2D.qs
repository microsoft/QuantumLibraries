// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

// FIXME: WIP for Ising All to all coupling model 
// ONe aspect at a time in each example. Keep it focused on one message. e.g. Ising for different coupling. Ising for different annealing times. Ising for different representation
namespace Microsoft.Quantum.Canon {
    open Microsoft.Quantum.Primitive
    

    /// # Summary
    /// Random real coupling with closed boundary conditions
    /// Use this to generate an array of couplings 
    operation GenerateRandom1DJCoupling(nSites: Int) : Double[] {
        body {
            mutable jCCouplings = new Double[nSites]
            let bitsOfRandomness = 8
            for(idx in 0..nSites-1) {
                set jCCouplings[idx] = RandomReal(bitsOfRandomness)
            }
            return jCCouplings
        }
    }
    function jC1DFromArray(data: Double[], schedule: Double, idx: Int): Double {
        return schedule * data[idx]
    }

   

    /// Ising all-to-all model
    function IsingA2AGeneratorIndex(nSites : Int, hx : Double,  jC: ((Int, Int) -> Double), idxHamiltonian : Int) : GeneratorIndex
    {
        // when idxHamiltonian is in [0, nSites - 1], apply transverse field "hx"
        // when idxHamiltonian is in [nSites, 2 * nSites - 1], apply Ising couplings jC(idxSite)
        mutable idxPauliString = new Int[0]
        mutable idxQubits = new Int[0]
        mutable coeff = ToDouble(0)

        if(idxHamiltonian < nSites){
            set idxPauliString = [1]
            set idxQubits = [idxHamiltonian]
            set coeff = - 1.0 * hx
        }
        elif(idxHamiltonian < nSites + nSites * nSites){
            let idx = idxHamiltonian - nSites
            let idxI = idx % nSites
            let idxJ = idx / nSites
            if(idxI < idxJ){
                set coeff = - 1.0 * jC(idxI, idxJ)
                set idxPauliString = [3; 3]
                set idxQubits = [idxI; idxJ]
            }
        }
        let generatorIndex = GeneratorIndex((idxPauliString,[coeff]),idxQubits)
        return generatorIndex
    }

    // This provides a full description of the Hamiltonian -- a map from a 1D index to terms in the Hamiltonian, and then to unitaries
    function IsingA2AEvolutionGenerator(nSites : Int, hx : Double,  jC: ((Int, Int) -> Double)) : EvolutionGenerator
    {
        let nTerms = nSites + nSites * nSites
        let generatorIndexFunction = IsingA2AGeneratorIndex(nSites, hx,  jC, _)
        let generatorSystem = GeneratorSystem(nTerms, generatorIndexFunction)
        let evolutionSet = PauliEvolutionSet()
        return EvolutionGenerator(evolutionSet, generatorSystem)
    }

    // This implements a trotter step that approximates the time-evolution operator
    function IsingA2ATrotterStepB(  nSites : Int, 
                                    hxSchedule : (Double -> Double), 
                                    jCSchedule: ((Double, Int, Int) -> Double), 
                                    schedule: Double, 
                                    trotterOrder: Int, 
                                    trotterStepSize: Double) : (Qubit[] => () : Adjoint, Controlled)
    {
        let hx = hxSchedule(schedule)
        let jC = jCSchedule(schedule, _, _)
        let evolutionGenerator = IsingA2AEvolutionGenerator(nSites, hx,  jC)
        return TrotterStep(evolutionGenerator, trotterOrder, trotterStepSize)
    }

    // This implements a trotter step for time-dependent evolution
    function IsingA2ASchedule(nSites : Int, hxSchedule : (Double -> Double), jCSchedule: ((Double, Int, Int) -> Double), trotterOrder: Int): EvolutionSchedule
    {
        let trotterStep = IsingA2ATrotterStepB(nSites, hxSchedule, jCSchedule, _, trotterOrder, _)
        return EvolutionSchedule(trotterStep)
    } 

    // Generate random couplings for 
    operation jC2DGenerateRandom(nSites: Int) : Double[][] {
        body {
            mutable jCCouplings = new Double[][nSites]
            let bitsOfRandomness = 8
            for(idxI in 0..nSites-1){
                mutable jCCouplingsJ = new Double[nSites]
                for(idxJ in 0..nSites-1){
                    set jCCouplingsJ[idxJ] = RandomReal(bitsOfRandomness)
                }
                set jCCouplings[idxI] = jCCouplingsJ
            }
            return jCCouplings
        }
    }
    function jC2DFromArray(data: Double[][], schedule: Double, idxI: Int, idxJ : Int): Double {
        return schedule * data[idxI][idxJ]
    }


    /// Adiabatic state preparation followed by phase estimation and measurement of sites
    /// This initializes the qubits in an easy-to-prepare eigenstate of the initial Hamiltonian
    operation IsingA2AStatePrep(qubits : Qubit[]) : (){
        body{
            ApplyToEachAC(X, qubits)
            ApplyToEachAC(H, qubits)
        }
        adjoint auto
        controlled auto
        controlled adjoint auto
    }

    //An example os the Ising model with random all-to-all couplings

    operation IsingA2ACriticalEstimate(nSites : Int, adiabaticTime: Double, bitsPrecision : Int, trotterOrder : Int, scheduleSteps : Int, qpeStepSize:Double) : (Double, Result[])
    {
        body{
            let hx =hxLinear( ToDouble(1), _)
            let jCData = jC2DGenerateRandom(nSites)
            let jC = jC2DFromArray(jCData, _, _, _)

            
            let statePrepUnitary = IsingA2AStatePrep
            let evolutionSchedule = IsingA2ASchedule(nSites, hx, jC, trotterOrder)
            ///let qpeUnitary = (evolutionSchedule)(ToDouble(1), qpeStepSize,_)
            let qpeUnitary = NoOp
            mutable phaseEst = ToDouble(0)
            mutable results = new Result[nSites]
            using(qubits = Qubit[nSites]){
                ResetAll(qubits)
                (AdiabaticStatePrep(adiabaticTime, scheduleSteps, evolutionSchedule, statePrepUnitary))(qubits)
                set phaseEst = RobustPhaseEstimation(bitsPrecision, OracleToDiscrete(qpeUnitary), qubits) / qpeStepSize  
                set results = MultiM(qubits)
                ResetAll(qubits)
            }

            return (phaseEst, results)
        }
    }

}

