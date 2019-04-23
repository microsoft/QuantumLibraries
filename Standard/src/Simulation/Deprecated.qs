// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {
    open Microsoft.Quantum.Simulation;
    open Microsoft.Quantum.Oracles;

    /// # Deprecated
    /// Please use @"microsoft.quantum.simuilation.pnormalized".
    operation AdiabaticStateEnergyUnitary(nQubits : Int, statePrepUnitary : (Qubit[] => Unit), adiabaticUnitary : (Qubit[] => Unit), qpeUnitary : (Qubit[] => Unit is Adj + Ctl), phaseEstAlgorithm : ((DiscreteOracle, Qubit[]) => Double)) : Double {
        _Renamed("Microsoft.Quantum.Canon.AdiabaticStateEnergyUnitary", "Microsoft.Quantum.Simulation.EstimateEnergyWithAdiabaticEvolution");
        return EstimateEnergyWithAdiabaticEvolution(nQubits, statePrepUnitary, adiabaticUnitary, qpeUnitary, phaseEstAlgorithm);
    }
}
