// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

// !! PLACEHOLDER METHODS !!
// Waiting until Microsoft.Quantum.Logical is PR'd
namespace Microsoft.Quantum.Diagnostics {
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Canon;

    /// # Summary
    /// Checks whether two ints are equal
    ///
    /// # Input
    /// ## first
    /// First int to check
    ///
    /// ## second
    /// Second int to check
    /// 
    /// # Output
    /// True if the ints are equal, False if the ints differ.
    function ClaimEqualInt(first : Int, second : Int) : Bool {
        if (first == second) {
            return true;
        }
        return false;
    }

    function ClaimDifferentInt(first : Int, second : Int) : Bool {
        if (first != second) {
            return true;
        }
        return false;
    }

}
