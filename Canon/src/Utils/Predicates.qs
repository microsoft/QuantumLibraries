// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon
{
    
    /// # Summary
    /// Tests if a given Result value is equal to `Zero`.
    ///
    /// # Input
    /// ## input
    /// `Result` value to be tested.
    /// # Output
    /// Returns `true` if `input` is equal to `Zero`.
    function IsResultZero (input : Result) : Bool
    {
        return input == Zero;
    }
    
    
    /// # Summary
    /// Tests if a given Result value is equal to `One`.
    ///
    /// # Input
    /// ## input
    /// `Result` value to be tested.
    /// # Output
    /// Returns `true` if `input` is equal to `One`.
    function IsResultOne (input : Result) : Bool
    {
        return input == One;
    }
    
}


