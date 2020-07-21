// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Diagnostics {
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Canon;


    operation AllowAtMostNCallsCA<'TInput, 'TOutput>(
        nTimes : Int, op : ('TInput => 'TOutput is Adj + Ctl)
    )
    : Unit is Adj {
    }

}
