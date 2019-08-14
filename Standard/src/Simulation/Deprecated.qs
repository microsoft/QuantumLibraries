// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {
    open Microsoft.Quantum.Simulation;
    open Microsoft.Quantum.Oracles;
    open Microsoft.Quantum.Warnings;
    open Microsoft.Quantum.Convert;

    /// # Deprecated
    /// Please use @"microsoft.quantum.simulation.estimateenergywithadiabaticevolution".
    operation AdiabaticStateEnergyUnitary(nQubits : Int, statePrepUnitary : (Qubit[] => Unit), adiabaticUnitary : (Qubit[] => Unit), qpeUnitary : (Qubit[] => Unit is Adj + Ctl), phaseEstAlgorithm : ((DiscreteOracle, Qubit[]) => Double)) : Double {
        _Renamed("Microsoft.Quantum.Canon.AdiabaticStateEnergyUnitary", "Microsoft.Quantum.Simulation.EstimateEnergyWithAdiabaticEvolution");
        return EstimateEnergyWithAdiabaticEvolution(nQubits, statePrepUnitary, adiabaticUnitary, qpeUnitary, phaseEstAlgorithm);
    }

    /// # Deprecated
    /// Please use @"microsoft.quantum.convert.intaspauli".
    function IntToPauli(idx : Int) : Pauli {
        _Renamed("Microsoft.Quantum.Simulation.IntToPauli", "Microsoft.Quantum.Convert.IntAsPauli");
        return IntAsPauli(idx);
    }

    /// # Deprecated
    /// Please use @"microsoft.quantum.convert.intarrayaspauliarray".
    function IntsToPaulis(ints : Int[]) : Pauli[] {
        _Renamed("Microsoft.Quantum.Simulation.IntsToPaulis", "Microsoft.Quantum.Convert.IntArrayAsPauliArray");
        return IntArrayAsPauliArray(ints);
    }

    /// # Deprecated
    /// Use `::NTerms` instead.
    function GetGeneratorSystemNTerms (generatorSystem : GeneratorSystem) : Int {
        _Removed(
            "Microsoft.Quantum.Simulation.GetGeneratorSystemNTerms",
            "Use ::NTerms instead."
        );
        return generatorSystem::NTerms;
    }


    /// # Deprecated
    /// Use `::Term` instead.
    function GetGeneratorSystemFunction (generatorSystem : GeneratorSystem) : (Int -> GeneratorIndex) {
        _Removed(
            "Microsoft.Quantum.Simulation.GetGeneratorSystemFunction",
            "Use ::Term instead."
        );
        return generatorSystem::Term;
    }


}
