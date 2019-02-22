// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {
    open Microsoft.Quantum.Primitive;

    /// # Summary
    /// When called from a deprecated function or operation, logs a message
    /// warning about the deprecation and informing the user of the new name.
    function Renamed(oldName : String, newName : String) : Unit {
        Message($"[WARNING] The callable {oldName} has been deprecated in favor of {newName}.");
    }

}
