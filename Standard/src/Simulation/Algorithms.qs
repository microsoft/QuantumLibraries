// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Simulation {
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Extensions.Convert;

    // A simulation technique converts an EvolutionGenerator to time evolution
    // by the encoded system for some time step
    // Here is an example of a simulation technique.
    
    /// # Summary
    /// Implements time-evolution by a term contained in a `GeneratorSystem`.
    ///
    /// # Input
    /// ## evolutionGenerator
    /// A complete description of the system to be simulated.
    /// ## idx
    /// Integer index to a term in the described system.
    /// ## stepsize
    /// Multiplier on duration of time-evolution by term indexed by `idx`.
    /// ## qubits
    /// Qubits acted on by simulation.
    operation TrotterStepImpl (evolutionGenerator : EvolutionGenerator, idx : Int, stepsize : Double, qubits : Qubit[]) : Unit is Adj + Ctl {
        let (evolutionSet, generatorSystem) = evolutionGenerator!;
        let (nTerms, generatorSystemFunction) = generatorSystem!;
        let generatorIndex = generatorSystemFunction(idx);
        (evolutionSet!(generatorIndex))!(stepsize, qubits);
    }
    
    
    /// # Summary
    /// Implements a single time-step of time-evolution by the system
    /// described in an `EvolutionGenerator` using a Trotter–Suzuki
    /// decomposition.
    ///
    /// # Input
    /// ## evolutionGenerator
    /// A complete description of the system to be simulated.
    /// ## trotterOrder
    /// Order of Trotter integrator. This must be either 1 or an even number.
    /// ## trotterStepSize
    /// Duration of simulated time-evolution in single Trotter step.
    ///
    /// # Output
    /// Unitary operation that approximates a single step of time-evolution
    /// for duration `trotterStepSize`.
    ///
    /// # Remarks
    /// For more on the Trotter–Suzuki decomposition, see
    /// [Time-Ordered Composition](/quantum/libraries/control-flow#time-ordered-composition).
    function TrotterStep (evolutionGenerator : EvolutionGenerator, trotterOrder : Int, trotterStepSize : Double) : (Qubit[] => Unit : Adjoint, Controlled)
    {
        let (evolutionSet, generatorSystem) = evolutionGenerator!;
        let (nTerms, generatorSystemFunction) = generatorSystem!;
        
        // The input to DecomposeIntoTimeStepsCA has signature
        // (Int, ((Int, Double, Qubit[]) => () : Adjoint, Controlled))
        let trotterForm = (nTerms, TrotterStepImpl(evolutionGenerator, _, _, _));
        return (DecomposeIntoTimeStepsCA(trotterForm, trotterOrder))(trotterStepSize, _);
    }
    
    
    // This simulation algorithm takes (timeMax, EvolutionGenerator,
    // register) and other algorithm-specific parameters (trotterStepSize,
    // trotterOrder), and performs evolution under the EvolutionGenerator
    // for time = timeMax.
    
    /// # Summary
    /// Makes repeated calls to `TrotterStep` to approximate the
    /// time-evolution operator exp(_-iHt_).
    ///
    /// # Input
    /// ## trotterStepSize
    /// Duration of simulated time-evolution in single Trotter step.
    /// ## trotterOrder
    /// Order of Trotter integrator. This must be either 1 or an even number.
    /// ## maxTime
    /// Total duration of simulation $t$.
    /// ## evolutionGenerator
    /// A complete description of the system to be simulated.
    /// ## qubits
    /// Qubits acted on by simulation.
    operation TrotterSimulationAlgorithmImpl (trotterStepSize : Double, trotterOrder : Int, maxTime : Double, evolutionGenerator : EvolutionGenerator, qubits : Qubit[]) : Unit
    {
        body (...)
        {
            let nTimeSlices = Microsoft.Quantum.Extensions.Math.Ceiling(maxTime / trotterStepSize);
            let resizedTrotterStepSize = maxTime / ToDouble(nTimeSlices);
            
            for (idxTimeSlice in 0 .. nTimeSlices - 1)
            {
                (TrotterStep(evolutionGenerator, trotterOrder, resizedTrotterStepSize))(qubits);
            }
        }
        
        adjoint invert;
        controlled distribute;
        controlled adjoint distribute;
    }
    
    
    /// # Summary
    /// `SimulationAlgorithm` function that uses a Trotter–Suzuki
    /// decomposition to approximate the time-evolution operator _exp(-iHt)_.
    ///
    /// # Input
    /// ## trotterStepSize
    /// Duration of simulated time-evolution in single Trotter step.
    /// ## trotterOrder
    /// Order of Trotter integrator. This must be either 1 or an even number.
    ///
    /// # Output
    /// A `SimulationAlgorithm` type.
    function TrotterSimulationAlgorithm (trotterStepSize : Double, trotterOrder : Int) : SimulationAlgorithm
    {
        return SimulationAlgorithm(TrotterSimulationAlgorithmImpl(trotterStepSize, trotterOrder, _, _, _));
    }
    
    
    // This simple time-dependent simulation algorithm implements a
    // sequence of uniformly-sized trotter steps
    
    /// # Summary
    /// Implementation of multiple Trotter steps to approximate a unitary
    /// operator that solves the time-dependent Schrödinger equation.
    ///
    /// # Input
    /// ## trotterStepSize
    /// Duration of simulated time-evolution in single Trotter step.
    /// ## trotterOrder
    /// Order of Trotter integrator. This must be either 1 or an even number.
    /// ## maxTime
    /// Total duration of simulation $t$.
    /// ## evolutionSchedule
    /// A complete description of the time-dependent system to be simulated.
    /// ## qubits
    /// Qubits acted on by simulation.
    operation TimeDependentTrotterSimulationAlgorithmImpl (trotterStepSize : Double, trotterOrder : Int, maxTime : Double, evolutionSchedule : EvolutionSchedule, qubits : Qubit[]) : Unit
    {
        body (...)
        {
            let nTimeSlices = Microsoft.Quantum.Extensions.Math.Ceiling(maxTime / trotterStepSize);
            let resizedTrotterStepSize = maxTime / ToDouble(nTimeSlices);
            
            for (idxTimeSlice in 0 .. nTimeSlices - 1)
            {
                let schedule = ToDouble(idxTimeSlice) / ToDouble(nTimeSlices);
                let (evolutionSet, generatorSystemTimeDependent) = evolutionSchedule!;
                let generatorSystem = generatorSystemTimeDependent(schedule);
                let evolutionGenerator = EvolutionGenerator(evolutionSet, generatorSystem);
                (TrotterSimulationAlgorithm(resizedTrotterStepSize, trotterOrder))!(resizedTrotterStepSize, evolutionGenerator, qubits);
            }
        }
        
        adjoint invert;
        controlled distribute;
        controlled adjoint distribute;
    }
    
    
    /// # Summary
    /// `TimeDependentSimulationAlgorithm` function that uses a Trotter–Suzuki
    /// decomposition to approximate a unitary operator that solves the
    /// time-dependent Schrodinger equation.
    ///
    /// # Input
    /// ## trotterStepSize
    /// Duration of simulated time-evolution in single Trotter step.
    /// ## trotterOrder
    /// Order of Trotter integrator. This must be either 1 or an even number.
    ///
    /// # Output
    /// A `TimeDependentSimulationAlgorithm` type.
    function TimeDependentTrotterSimulationAlgorithm (trotterStepSize : Double, trotterOrder : Int) : TimeDependentSimulationAlgorithm
    {
        return TimeDependentSimulationAlgorithm(TimeDependentTrotterSimulationAlgorithmImpl(trotterStepSize, trotterOrder, _, _, _));
    }
    
}


