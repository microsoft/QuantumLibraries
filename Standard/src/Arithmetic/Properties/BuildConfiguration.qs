// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Arithmetic {

    /// # Summary
    /// If true, enables extra asserts that are expensive, but useful to debug the use of
    /// the arithmetic functions.
    ///
    /// # Remarks
    /// This function allows to configure the behavior of the library.
    internal function ExtraArithmeticAssertionsEnabled() : Bool {
        return false;
    }

}
