// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {

    function Fst<'T, 'U>(pair : ('T, 'U)) : 'T {
        let (fst, snd) = pair
        return fst
    }

    function Snd<'T, 'U>(pair : ('T, 'U)) : 'U {
        let (fst, snd) = pair
        return snd
    }

}
