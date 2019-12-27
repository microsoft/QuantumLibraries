// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.AmplitudeAmplification {
    open Microsoft.Quantum.Convert;
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Math;

    /// # Summary
    /// Computes partial reflection phases for standard amplitude
    /// amplification.
    ///
    /// # Input
    /// ## nIterations
    /// Number of amplitude amplification iterations to generate partial
    /// reflection phases for.
    ///
    /// # Output
    /// An operation that implements phases specified as partial reflections
    ///
    /// # Remarks
    /// All phases are $\pi$, except for the first reflection about the start
    /// state and the last reflection about the target state, which are $0$.
    function StandardReflectionPhases(nIterations : Int) : ReflectionPhases {
        let commonPhases = ConstantArray(nIterations, PI());
        let targetPhases = commonPhases + [0.0];
        let startPhases = [0.0] + commonPhases;
        return ReflectionPhases(startPhases, targetPhases);
    }

    // We use the phases in "Fixed-Point Amplitude Amplification with an
    // Optimal Number of Queires" [YoderLowChuang2014]
    // See also "Methodology of composite quantum gates" [LowYoderChuang2016]
    // for phases in the `RotationPhases` format

    /// # Summary
    /// Computes partial reflection phases for fixed-point amplitude
    /// amplification.
    ///
    /// # Input
    /// ## nQueries
    /// Number of queries to the state preparation oracle. Must be an odd
    /// integer.
    /// ## successMin
    /// Target minimum success probability.
    ///
    /// # Output
    /// Array of phases that can be used in fixed-point amplitude amplification
    /// quantum algorithm implementation.
    ///
    /// # References
    /// We use the phases in "Fixed-Point Amplitude Amplification with
    /// an Optimal Number of Queries"
    /// - [YoderLowChuang2014](https://arxiv.org/abs/1409.3305)
    /// See also "Methodology of composite quantum gates"
    /// - [LowYoderChuang2016](https://arxiv.org/abs/1603.03996)
    /// for phases in the `RotationPhases` format.
    function FixedPointReflectionPhases(nQueries : Int, successMin : Double)
    : ReflectionPhases {
        let twoPi = 2.0 * PI();
        mutable phasesRot = new Double[nQueries];
        let nQueriesDouble = IntAsDouble(nQueries);
        set phasesRot w/= 0 <- 0.0;
        let beta = Cosh((1.0 / nQueriesDouble) * ArcCosh(Sqrt(successMin)));
        let alpha = Sqrt(1.0 - beta * beta);

        for (idxPhases in 1 .. nQueries - 1) {
            set phasesRot w/= idxPhases <-
                phasesRot[idxPhases - 1] +
                2.0 * ArcTan(
                    Tan(twoPi * IntAsDouble(idxPhases) / nQueriesDouble) * alpha
                );
        }

        return RotationPhasesAsReflectionPhases(RotationPhases(phasesRot));
    }

}
