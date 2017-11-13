// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Tests {
    open Microsoft.Quantum.Canon;    

    function ResultAsIntTest() : () {
        // FIXME: Write new Assert for ints, or better yet, use generics.
        if (ResultAsInt([Zero; Zero]) != 0) {
            fail "Expected 0.";
        }
        if (ResultAsInt([One; Zero]) != 1) {
            fail "Expected 1.";
        }
        if (ResultAsInt([Zero; One]) != 2) {
            fail "Expected 2.";
        }
        if (ResultAsInt([One; One]) != 3) {
            fail "Expected 3.";
        }
    }
}
