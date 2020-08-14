// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Math {
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Convert;
    open Microsoft.Quantum.Random;

    @Deprecated("Microsoft.Quantum.Random.DrawRandomInt")
    operation RandomIntPow2 (maxBits : Int) : Int {
        return DrawRandomInt(0, 2^maxBits - 1);
    }
    
    @Deprecated("Microsoft.Quantum.Random.DrawRandomInt")
    operation RandomInt (maxInt : Int) : Int {
        return DrawRandomInt(0, maxInt - 1);
    }
    
    
    @Deprecated("Microsoft.Quantum.Random.DrawRandomDouble")
    operation RandomReal (bitsRandom : Int) : Double {
        if (bitsRandom < 1) {
            fail $"Number of random bits must be greater than 0.";
        }
        
        return IntAsDouble(RandomIntPow2(bitsRandom)) / PowD(2.0, IntAsDouble(bitsRandom));
    }

    @Deprecated("Microsoft.Quantum.Random.DrawRandomPauli")
    operation RandomSingleQubitPauli() : Pauli {
        let probs = [0.5, 0.5, 0.5, 0.5];
        let idxPauli = Random(probs);
        let singleQubitPaulis = [PauliI, PauliX, PauliY, PauliZ];
        return singleQubitPaulis[idxPauli];
    }

}


