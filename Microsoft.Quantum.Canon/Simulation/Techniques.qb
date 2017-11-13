// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {
    open Microsoft.Quantum.Extensions.Convert;
    open Microsoft.Quantum.Extensions.Math;

    // Private operation interpolating between two generators.
    // Not documented independently of AdiabaticEvolution, below.
    operation _AdiabaticEvolution(
            adiabaticTime: Double, 
            evolutionGeneratorStart: EvolutionGenerator,
            evolutionGeneratorEnd: EvolutionGenerator,
            simulationAlgorithmTimeDependent: SimulationAlgorithmTimeDependent,
            qubits: Qubit[]) : () 
    {
        body {
            // evolutionSetStart and evolutionSetEnd must be identical
            let (evolutionSetStart, generatorSystemStart) = evolutionGeneratorStart;
            let (evolutionSetEnd, generatorSystemEnd) = evolutionGeneratorEnd;
            let generatorSystemTimeDependent = InterpolateGeneratorSystems(generatorSystemStart, generatorSystemEnd);
            let evolutionSchedule = EvolutionSchedule(evolutionSetStart, generatorSystemTimeDependent);
            simulationAlgorithmTimeDependent(adiabaticTime, evolutionSchedule, qubits);
        }
        adjoint auto
        controlled auto
        controlled adjoint auto
    }

    //FIXME EXEC : error QB5001: internal error in code gen for partial application
    //FIXME AdiabaticEvolution should be implemented in this manner without requiring _AdiabaticEvolution
    //function FIXMEAdiabaticEvolution(   adiabaticTime: Double, 
    //                                   evolutionGeneratorStart: EvolutionGenerator,
    //                                   evolutionGeneratorEnd: EvolutionGenerator,
    //                                   simulationAlgorithmTimeDependent: SimulationAlgorithmTimeDependent) : (Qubit[] => () : Adjoint, Controlled) {
    //
        // evolutionSetStart and evolutionSetEnd must be identical
    //    let (evolutionSetStart, generatorSystemStart) = evolutionGeneratorStart
    //    let (evolutionSetEnd, generatorSystemEnd) = evolutionGeneratorEnd
    //    let generatorSystemTimeDependent = InterpolateGeneratorSystems(generatorSystemStart, generatorSystemEnd)
    //    let evolutionSchedule = EvolutionSchedule(evolutionSetStart, generatorSystemTimeDependent)
    //    return simulationAlgorithmTimeDependent(adiabaticTime, evolutionSchedule, _)
    //}

    // FIXME: (The following is CG's opinion, please feel free to disregard accordingly.)
    //        This should be renamed, as it's not in general adiabatic.
    //        In particular, someone could pass arbitrary evolution times,
    //        rather than one meeting the adiabtic conditions.
    //        Similarly, arguments to this function should be renamed, I think.
    /// # Summary
    /// This interpolates between two generators with a uniform schedule,
    /// returning an operation that applies simulated evolution under
    /// the resulting time-dependent generator to a qubit register.
    ///
    /// # Input
    /// ## adiabaticTime
    /// Time to perform the interpolation over.
    /// ## evolutionGeneratorStart
    /// Initial generator to simulate evolution under.
    /// ## evolutionGeneratorEnd
    /// Final generator to simulate evolution under.
    /// ## simulationAlgorithmTimeDependent
    /// A time-dependent simulation algorithm that will be used
    /// to simulate evolution during the uniform interpolation schedule.
    ///
    /// # Remarks
    /// If the interpolation time is chosen to meet the adiabatic conditions,
    /// then this function returns an operation which performs adiabatic
    /// state preparation for the final dynamical generator.
    ///
    /// # References
    /// **TODO** This would be a good place to link to Nathan's paper on
    /// adiabatic conditions.
    function AdiabaticEvolution(
        adiabaticTime: Double, 
        evolutionGeneratorStart: EvolutionGenerator,
        evolutionGeneratorEnd: EvolutionGenerator,
        simulationAlgorithmTimeDependent: SimulationAlgorithmTimeDependent)
        : (Qubit[] => () : Adjoint, Controlled)
    {  
        return _AdiabaticEvolution(
            adiabaticTime, evolutionGeneratorStart, evolutionGeneratorEnd,
            simulationAlgorithmTimeDependent, _
        );
    }



    /// summary:
    ///     Performs state preparation, and adiabatic state preparation, and phase est.
    operation AdiabaticStateEnergyUnitary(  statePrepUnitary: (Qubit[] => ()),
                                            adiabaticUnitary: (Qubit[] => ()),
                                            qpeUnitary: (Qubit[] => () :  Adjoint, Controlled),
                                            phaseEstAlgorithm : ((DiscreteOracle, Qubit[]) => Double),
                                            qubits: Qubit[]) : Double {
        body {
            statePrepUnitary(qubits);
            adiabaticUnitary(qubits);
            let phaseEst = phaseEstAlgorithm(OracleToDiscrete(qpeUnitary), qubits);
            return phaseEst;
        }
    }

    //Qubits + State prep + Adiabatic state prep + energy estimation + measure in basis
    operation AdiabaticStateEnergyEstimate( nQubits : Int, 
                                            statePrepUnitary: (Qubit[] => ()),
                                            adiabaticUnitary: (Qubit[] => ()),
                                            qpeUnitary: (Qubit[] => () :  Adjoint, Controlled),
                                            phaseEstAlgorithm : ((DiscreteOracle, Qubit[]) => Double)) : Double {
        body {
            mutable phaseEst = ToDouble(0);
            using (qubits = Qubit[nQubits]) {
                set phaseEst = AdiabaticStateEnergyUnitary( statePrepUnitary, adiabaticUnitary, qpeUnitary, phaseEstAlgorithm, qubits);
                ResetAll(qubits);
            }
            return phaseEst;

        }
    }

    //Qubits + State Prep +  quantum phase estimation
    operation EstimateEnergy(nQubits : Int,
                             statePrepUnitary: (Qubit[] => () ),
                             qpeUnitary: (Qubit[] => () :  Adjoint, Controlled),
                             phaseEstAlgorithm : ((DiscreteOracle, Qubit[]) => Double) ) : Double {
         body {
            let phaseEst = AdiabaticStateEnergyEstimate( nQubits, statePrepUnitary, NoOp, qpeUnitary, phaseEstAlgorithm);
            return phaseEst;
        }
    }



//For LCU
}
