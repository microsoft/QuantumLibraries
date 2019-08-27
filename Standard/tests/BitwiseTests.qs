// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Tests {
    open Microsoft.Quantum.Bitwise;
    open Microsoft.Quantum.Diagnostics;

    function ShiftTest() : Unit {
        let smallValue = 5; // 0b101
        EqualityFactI(20, LeftShiftedI(smallValue, 2), "Shifted values incorrect.");
        EqualityFactI(4, RightShiftedI(smallValue, 1), "Shifted values incorrect.");

        // 0b‭1111111011011100101110101001100001110110010101000011001000010000‬
        let bigValue = 0xdb97530eca864200L;
        EqualityFactL(
            0xffdb97530eca8642L,
            LeftShiftedL(bigValue, 5),
            "Shfited values incorrect."
        );
        EqualityFactL(
            0xffdb97530eca8642L,
            RightShiftedL(bigValue, 3),
            "Shfited values incorrect."
        );
    }

}
