// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.AmplitudeAmplification {
    open Microsoft.Quantum.Oracles;


    /// # Deprecated
    /// Please use
    /// @"microsoft.quantum.amplitudeamplification.rotationphasesasreflectionphases".
    @Deprecated("Microsoft.Quantum.AmplitudeAmplification.RotationPhasesAsReflectionPhases")
    function AmpAmpRotationToReflectionPhases (rotPhases : RotationPhases)
    : ReflectionPhases {
        return RotationPhasesAsReflectionPhases(rotPhases);
    }

    /// # Deprecated
    /// Please use
    /// @"microsoft.quantum.amplitudeamplification.standardreflectionphases".
    @Deprecated("Microsoft.Quantum.AmplitudeAmplification.StandardReflectionPhases")
    function AmpAmpPhasesStandard(nIterations : Int) : ReflectionPhases {
        return StandardReflectionPhases(nIterations);
    }

    /// # Deprecated
    /// Please use
    /// @"microsoft.quantum.amplitudeamplification.fixedpointreflectionphases".
    @Deprecated("Microsoft.Quantum.AmplitudeAmplification.FixedPointReflectionPhases")
    function AmpAmpPhasesFixedPoint(nQueries : Int, successMin : Double) : ReflectionPhases {
        return FixedPointReflectionPhases(nQueries, successMin);
    }

    /// # Deprecated
    /// Please use @"microsoft.quantum.amplitudeamplification.obliviousamplitudeamplificationfrompartialreflections".
    @Deprecated("Microsoft.Quantum.AmplitudeAmplification.ObliviousAmplitudeAmplificationFromPartialReflections")
    function AmpAmpObliviousByReflectionPhases(
        phases : ReflectionPhases,
        startStateReflection : ReflectionOracle,
        targetStateReflection : ReflectionOracle,
        signalOracle : ObliviousOracle
    )
    : ((Qubit[], Qubit[]) => Unit is Adj + Ctl) {
        return ObliviousAmplitudeAmplificationFromPartialReflections(
            phases, startStateReflection, targetStateReflection, signalOracle
        );
    }

    /// # Deprecated
    /// Please use @"microsoft.quantum.amplitudeamplification.obliviousamplitudeamplificationfromstatepreparation".
    @Deprecated("Microsoft.Quantum.AmplitudeAmplification.ObliviousAmplitudeAmplificationFromStatePreparation")
    function AmpAmpObliviousByOraclePhases(
        phases : ReflectionPhases,
        startStateOracle : DeterministicStateOracle,
        signalOracle : ObliviousOracle,
        idxFlagQubit : Int
    )
    : ((Qubit[], Qubit[]) => Unit is Adj + Ctl) {
        return ObliviousAmplitudeAmplificationFromStatePreparation(
            phases, startStateOracle, signalOracle, idxFlagQubit
        );
    }

    /// # Deprecated
    /// Please use @"microsoft.quantum.amplitudeamplification.amplitudeamplificationfrompartialreflections".
    @Deprecated("Microsoft.Quantum.AmplitudeAmplification.AmplitudeAmplificationFromPartialReflections")
    function AmpAmpByReflectionPhases(
        phases : ReflectionPhases,
        startStateReflection : ReflectionOracle,
        targetStateReflection : ReflectionOracle
    )
    : (Qubit[] => Unit is Adj + Ctl) {
        return AmplitudeAmplificationFromPartialReflections(
            phases, startStateReflection, targetStateReflection
        );
    }

    /// # Deprecated
    /// Please use @"microsoft.quantum.amplitudeamplification.amplitudeamplificationfromstatepreparation".
    @Deprecated("Microsoft.Quantum.AmplitudeAmplification.AmplitudeAmplificationFromStatePreparation")
    function AmpAmpByOraclePhases(
        phases : ReflectionPhases,
        stateOracle : StateOracle,
        idxFlagQubit : Int
    )
    : (Qubit[] => Unit is Adj + Ctl) {
        return AmplitudeAmplificationFromStatePreparation(
            phases, stateOracle, idxFlagQubit
        );
    }

    /// # Deprecated
    /// Please use @"microsoft.quantum.amplitudeamplification.standardamplitudeamplification".
    @Deprecated("Microsoft.Quantum.AmplitudeAmplification.StandardAmplitudeAmplification")
    function AmpAmpByOracle(
        nIterations : Int,
        stateOracle : StateOracle,
        idxFlagQubit : Int
    )
    : (Qubit[] => Unit is Adj + Ctl) {
        return StandardAmplitudeAmplification(nIterations, stateOracle, idxFlagQubit);
    }


    /// # Deprecated
    /// Please use @"microsoft.quantum.amplitudeamplification.applyfixedpointamplification".
    @Deprecated("Microsoft.Quantum.AmplitudeAmplification.ApplyFixedPointAmplification")
    operation AmpAmpRUSByOracle(statePrepOracle : StateOracle, startQubits : Qubit[])
    : Unit {
        ApplyFixedPointAmplification(statePrepOracle, startQubits);
    }

}
