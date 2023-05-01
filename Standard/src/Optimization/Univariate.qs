// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Optimization {
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Math;

    /// # Summary
    /// Represents the result of optimizing a univariate function.
    ///
    /// # Named Items
    /// ## Coordinate
    /// Input at which an optimum was found.
    /// ## Value
    /// Value returned by the function at its optimum.
    /// ## NQueries
    /// The amount of times the function was called.
    newtype UnivariateOptimizationResult = (
        Coordinate : Double,
        Value : Double,
        NQueries : Int
    );

    /// # Summary
    /// Returns the width of an interval.
    internal function Width(left : Double, right : Double) : Double {
        Fact(left <= right, "Left endpoint of bounds must be less than or equal to right endpoint.");
        return right - left;
    }

    /// # Summary
    /// Given an interval, returns two probes that contract the given
    /// interval by a factor of the golden ratio.
    internal function NextProbes(left : Double, right : Double) : (Double, Double) {
        let goldenRatio = (Sqrt(5.0) + 1.0) / 2.0;
        let delta = (Width(left, right)) / goldenRatio;
        return (
            right - delta, left + delta
        );
    }

    /// # Summary
    /// Evaluates the given function at a coordinate, and returns the
    /// corresponding point.
    internal function ProbeValue(
        fn : (Double -> Double),
        coord : Double
    ) : (Double, Double) {
        return (coord, fn(coord));
    }

    /// # Summary
    /// Returns some local minimum for a univariate function over a bounded interval,
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
        Fact(tolerance > 0.0, "The tolerance value must be positive.");
        mutable interval = bounds;
        mutable leftProbe = ProbeValue(fn, Fst(NextProbes(interval)));
        mutable rightProbe = ProbeValue(fn, Snd(NextProbes(interval)));
        mutable queryAmount = 2;
        while (Width(interval) > tolerance) {
            if (Snd(leftProbe) < Snd(rightProbe)) {
                set interval = (Fst(interval), Fst(rightProbe));
                set rightProbe = leftProbe;
                set leftProbe = ProbeValue(fn, Fst(NextProbes(interval)));
            } else {
                set interval = (Fst(leftProbe), Snd(interval));
                set leftProbe = rightProbe;
                set rightProbe = ProbeValue(fn, Snd(NextProbes(interval)));
            }
            set queryAmount += 1;
        }

        // Return from the existing probes to avoid extra call to fn
        let result = Snd(leftProbe) < Snd(rightProbe) ? leftProbe | rightProbe;
        return UnivariateOptimizationResult(
            Fst(result), Snd(result), queryAmount
        );
    }

}
