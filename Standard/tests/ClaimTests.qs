// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Tests {
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Diagnostics;


    operation ClaimEqualIntTest() : Unit {
        let int1 = 10;
        let int2 = 24;

        EqualityFactB(ClaimEqualInt(int1, int1), true, "ClaimEqualInt says two ints are different when they are the same");
        EqualityFactB(ClaimEqualInt(int1, int2), false, "ClaimEqualInt says two ints are the same when they are different");
    }

}
