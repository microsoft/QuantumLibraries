// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {
    // Overview of simulation library.

    /// summary:
    ///     Represents a unitary operator.
    newtype Unitary = (Qubit[] => () : Adjoint, Controlled);

    /// summary:
    ///     Represents a unitary time-evolution operator. The first parameter is
    ///     is duration of time-evoltion, and the second parameter is the qubit 
    ///     register acted upon by the unitary
    newtype EvolutionUnitary = ((Double ,Qubit[]) => () : Adjoint, Controlled);

    /// summary:
    ///     Represents a set of gates that can be readily implemented and used
    ///     to implement simulation algorithms. Elements in the set are indexed 
    ///     by a "GeneratorIndex", and each set is described by a function
    ///     from "GeneratorIndex" to "EvolutionUnitary", which are operations 
    ///     parameterized by a real number representing time
    newtype EvolutionSet = (GeneratorIndex -> EvolutionUnitary);

    /// summary:
    ///     Represents a dynamical generator as a set of simulatable gates and
    ///     an expansion in terms of that basis.
    ///     Last parameter for number of terms
    newtype EvolutionGenerator = (EvolutionSet, GeneratorSystem);

    /// summary:
    ///     Represents a time-dependent dynamical generator. The Double 
    ///     parameter is a schedule in [0,1]
    newtype EvolutionSchedule = (EvolutionSet, (Double -> GeneratorSystem));

    /// summary:
    ///     A time-independent simulation technique converts an 
    ///     EvolutionGenerator to unitary time evolution for some time-interval
    newtype SimulationAlgorithm = ((Double, EvolutionGenerator, Qubit[]) => () : Adjoint, Controlled);

    /// summary:
    ///     A time-dependent simulation technique converts an EvolutionSchedule
    ///     to unitary time evolution for some time-interval
    newtype SimulationAlgorithmTimeDependent = ((Double, EvolutionSchedule, Qubit[]) => () : Adjoint, Controlled);

}
