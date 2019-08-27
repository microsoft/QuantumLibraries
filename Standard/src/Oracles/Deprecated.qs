// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Oracles {
    open Microsoft.Quantum.Warnings;

    /// # Deprecated
    /// Please use @"microsoft.quantum.amplitudeamplification.boundwithdeterministicstateoracle".
    function ObliviousOracleFromDeterministicStateOracle(ancillaOracle : DeterministicStateOracle, signalOracle : ObliviousOracle) : ObliviousOracle {
        _Renamed(
            "Microsoft.Quantum.Oracles.ObliviousOracleFromDeterministicStateOracle",
            "Microsoft.Quantum.Oracles.BoundWithDeterministicStateOracle"
        );
        return BoundWithDeterministicStateOracle(ancillaOracle, signalOracle);
    }

}

