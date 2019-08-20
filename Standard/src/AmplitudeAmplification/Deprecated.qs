// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.AmplitudeAmplification {
    open Microsoft.Quantum.Warnings;
    open Microsoft.Quantum.Math;

    /// # Deprecated
    /// Please use @"microsoft.quantum.amplitudeamplification.pnormalized".
    function AmpAmpRotationToReflectionPhases(rotPhases : RotationPhases) : ReflectionPhases {
        _Renamed(
            "Microsoft.Quantum.AmplitudeAmplification.AmpAmpRotationToReflectionPhases",
            "Microsoft.Quantum.AmplitudeAmplification.RotationPhasesAsReflectionPhases"
        );
        return RotationPhasesAsReflectionPhases(rotPhases);
    }

}

