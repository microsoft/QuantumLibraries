// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Simulation {

    // For an overview of the simulation library, see [Hamiltonian
    // Simulation](applications#hamiltonian-simulation)

    /// # Summary
	/// Represents a time-independent simulation algorithm.
	/// 
    /// A time-independent simulation technique converts an
    /// <xref:Microsoft.Quantum.Simulation.EvolutionGenerator>
    /// to unitary time evolution for some time-interval.
    ///
    /// # Description
    /// The inputs into the callable are:
    /// - The time interval of simulation.
    /// - A representation of the generator of dynamic evolution.
    /// - A register encoding the state of the system.
    ///
    /// # Example
    /// To apply the Trotter–Suzuki simulation algorithm to a register of
    /// qubits:
    ///
    /// ```qsharp
    /// operation EvolveUnderGenerator(generator : EvolutionGenerator, time : Double, register : Qubit[])
    /// : Unit is Adj + Ctl {
    ///     let trotterStepSize = 0.1;
    ///     let trotterOrder = 1;
    ///     let evolveFor = (TrotterSimulationAlgorithm(trotterStepSize, trotterOrder))!;
    ///     evolveFor(time, generator, register);
    /// }
    /// ```
    ///
    /// # See Also
    /// - Microsoft.Quantum.Simulation.TimeDependentSimulationAlgorithm
    newtype SimulationAlgorithm = ((Double, EvolutionGenerator, Qubit[]) => Unit is Adj + Ctl);

    /// # Summary
	/// Represents a time-dependent simulation algorithm.
	/// 
    /// A time-dependent simulation technique converts an
    /// <xref:Microsoft.Quantum.Simulation.EvolutionSchedule>
    /// to unitary time-evolution for some time-interval.
    ///
    /// # Description
    /// The inputs into the callable are:
    /// - The time interval of simulation.
    /// - A schedule mapping evolution time to the generator at that time.
    /// - A register encoding the state of the system.
    ///
    /// # Example
    /// To apply the Trotter–Suzuki simulation algorithm to a register of
    /// qubits:
    ///
    /// ```qsharp
    /// operation EvolveUnderTimeDependentGenerator(schedule : EvolutionSchedule, time : Double, register : Qubit[])
    /// : Unit is Adj + Ctl {
    ///     let trotterStepSize = 0.1;
    ///     let trotterOrder = 1;
    ///     let evolveFor = (TimeDependentTrotterSimulationAlgorithm(trotterStepSize, trotterOrder))!;
    ///     evolveFor(time, schedule, register);
    /// }
    /// ```
    ///
    /// # See Also
    /// - Microsoft.Quantum.Simulation.TimeDependentSimulationAlgorithm
    ///
    /// # See Also
    /// - Microsoft.Quantum.Simulation.SimulationAlgorithm
    newtype TimeDependentSimulationAlgorithm = ((Double, EvolutionSchedule, Qubit[]) => Unit is Adj + Ctl);
    
}


