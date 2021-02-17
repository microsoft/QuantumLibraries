// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Random {
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Synthesis;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Canon;

    // TODO
    operation DrawRandomSingleQubitClifford() : SingleQubitClifford {
        return Identity1C()
            w/ E <- DrawRandomInt(0, 2)
            w/ S <- DrawRandomInt(0, 3)
            w/ X <- DrawRandomInt(0, 1)
            w/ Omega <- DrawRandomInt(0, 7);
    }
}
