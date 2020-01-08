// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {

    function _Compose<'T, 'U, 'V> (outer : ('U -> 'V), inner : ('T -> 'U), target : 'T) : 'V {
        return outer(inner(target));
    }

    /// # Summary
	/// Returns the composition of two functions.
    ///
    /// # Description
    /// Given two functions $f$ and $g$, returns a new function representing
    /// $f \circ g$.
    ///
    /// # Input
    /// ## outer
    /// The second function to be applied.
    /// ## inner
    /// The first function to be applied.
    ///
    /// # Output
    /// A new function $h$ such that for all inputs $x$, $f(g(x)) = h(x)$.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The input type of the first function to be applied.
    /// ## 'U
    /// The output type of the first function to be applied and the input type
    /// of the second function to be applied.
    /// ## 'V
    /// The output type of the second function to be applied.
    function Compose<'T, 'U, 'V> (outer : ('U -> 'V), inner : ('T -> 'U)) : ('T -> 'V) {
        return _Compose(outer, inner, _);
    }

}
