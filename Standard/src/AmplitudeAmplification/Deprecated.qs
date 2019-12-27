// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.AmplitudeAmplification {


    /// # Deprecated
    /// Please use
    // @"microsoft.quantum.amplitudeamplification.rotationphasesasreflectionphases".
    @Deprecated("Microsoft.Quantum.AmplitudeAmplification.RotationPhasesAsReflectionPhases")
    function AmpAmpRotationToReflectionPhases (rotPhases : RotationPhases)
    : ReflectionPhases {
        return RotationPhasesAsReflectionPhases(rotPhases);
    }

    /// # Deprecated
    /// Please use
    // @"microsoft.quantum.amplitudeamplification.standardreflectionphases".
    @Deprecated("Microsoft.Quantum.AmplitudeAmplification.StandardReflectionPhases")
    function AmpAmpPhasesStandard(nIterations : Int) : ReflectionPhases {
        return StandardReflectionPhases(nIterations);
    }

    /// # Deprecated
    /// Please use
    // @"microsoft.quantum.amplitudeamplification.fixedpointreflectionphases".
    @Deprecated("Microsoft.Quantum.AmplitudeAmplification.FixedPointReflectionPhases")
    function AmpAmpPhasesFixedPoint(nQueries : Int, successMin : Double) : ReflectionPhases {
        return FixedPointReflectionPhases(nQueries, successMin);
    }

}
