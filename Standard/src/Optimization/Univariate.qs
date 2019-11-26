// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Optimization {
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Math;

    /// # Summary
    /// Represents the result of optimizing a univariate function.
    ///
    /// # Input
    /// ## Coordinate
    /// Input at which an optimium was found.
    /// ## Value
    /// Value returned by the function at its optimum.
    newtype UnivariateOptimizationResult = (
        Coordinate : Double,
        Value : Double
    );

    /// # Summary
    /// Returns the width of an interval.
    function _Width(left : Double, right : Double) : Double {
        return right - left;
    }

    /// # Summary
    /// Given an interval, returns a probe interval that contracts the given
    /// interval by a factor of the golden ratio.
    function _Probe(left : Double, right : Double) : (Double, Double) {
        let goldenRatio = (Sqrt(5.0) + 1.0) / 2.0;
        let delta = (_Width(left, right)) / goldenRatio;
        return (
            left + delta, right - delta
        );
    }

    /// # Summary
    /// Returns the midpoint for an interval.
    function _Midpoint(left : Double, right : Double) : (Double) {
        return (left + right) / 2.0;
    }

    /// # Summary
    /// Returns the local minimum for a univariate function over a bounded interval,
    /// using a golden interval search.
    ///
    /// # Input
    /// ## fn
    /// The univariate function to be minimized.
    /// ## bounds
    /// The interval in which the local minimum is to be found.
    /// ## tolerance
    /// The accuracy in the function input to be tolerated.
    /// The search for local optima will continue until the interval is
    /// smaller in width than this tolerance.
    ///
    /// # Output
    /// A local optima of the given function, located within the given bounds.
    function LocalUnivariateMinimum(
        fn : (Double -> Double),
        bounds : (Double, Double),
        tolerance : Double
    ) : UnivariateOptimizationResult {
        
        mutable interval = bounds;
        mutable probe = _Probe(interval);
        
        while (_Width(probe) > tolerance) {
            set interval =
                fn(Fst(probe)) < fn(Snd(probe))
                ? (Fst(interval), Snd(probe))
                | (Fst(probe), Snd(interval));
            set probe = _Probe(interval);
        }

        let mid = _Midpoint(interval);
        return UnivariateOptimizationResult(
            mid, fn(mid)
        );

    }

}
