// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.AmplitudeAmplification {
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Convert;
    open Microsoft.Quantum.Diagnostics;
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
        // In this implementation `nQueries` corresponds to $L$ in
        // arXiv:1409.3305.
        Fact(nQueries % 2 == 1, "nQueries must be odd");

        // Initializes L rotation phases, this also initializes the first
        // rotation phase with 0.0.
        mutable phasesRot = [0.0, size = nQueries];
        let nQueriesDouble = IntAsDouble(nQueries);

        // The success probability `successMin` is $1 - \delta^2$ in
        // arXiv:1409.3305. Variable `beta` corresponds to $\gamma^{-1}$ in
        // arXiv:1409.3305, right below Eq. (11)
        let beta = Cosh((1.0 / nQueriesDouble) * ArcCosh(Sqrt(1.0 / (1.0 - successMin))));

        // `alpha` is $\sqrt(1 - \gamma^2)$ in Eq. (11) in arXiv:1409.3305,
        // therefore it is $\sqrt(1 - (1 / \beta^2))$
        let alpha = Sqrt(1.0 - 1.0 / (beta * beta));

        // Iterative computation of rotation phases is described in Eq. (30) in
        // arXiv:1603.03996.  In there, we can set $j = 1$.
        let twoPi = 2.0 * PI();
        for idxPhases in 1 .. nQueries - 1 {
            set phasesRot w/= idxPhases <-
                phasesRot[idxPhases - 1] +
                2.0 * ArcTan(
                    Tan(twoPi * IntAsDouble(idxPhases) / nQueriesDouble) * alpha
                );
        }

        return RotationPhasesAsReflectionPhases(RotationPhases(phasesRot));
    }

}
