// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {

    // For an overview of the simulation library, see [Hamiltonian 
    // Simulation](applications#hamiltonian-simulation)

    /// # Summary
    /// A time-independent simulation technique converts an
    //  @"microsoft.quantum.canon.evolutiongenerator"
    /// to unitary time evolution for some time-interval.
    newtype SimulationAlgorithm = ((Double, EvolutionGenerator, Qubit[]) => () : Adjoint, Controlled);

    /// # Summary
    /// A time-dependent simulation technique converts an
    /// @"microsoft.quantum.canon.evolutionschedule"
    /// to unitary time-evolution for some time-interval.
    newtype TimeDependentSimulationAlgorithm = ((Double, EvolutionSchedule, Qubit[]) => () : Adjoint, Controlled);

}
