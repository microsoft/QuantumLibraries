// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {
    /// # Summary
    /// The `IsResultZero` function tests if a given Result value is equal to `Zero`.     
    function IsResultZero (input : Result) : Bool {
        return (input == Zero);
    }

    /// # Summary
    /// The `IsResultOne` function tests if a given Result value is equal to `One`.
    function IsResultOne (input : Result) : Bool {
        return (input == One);
    }
    
}
