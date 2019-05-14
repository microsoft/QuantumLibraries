// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.
namespace Microsoft.Quantum.Tests {
    open Microsoft.Quantum.Math;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Canon;
    
    
    /// # Summary
    /// Checks that @"microsoft.quantum.canon.randomint" obeys ranges.
    operation RandomIntRangeTest () : Unit {
        let randomInt = RandomInt(45);
        if (randomInt > 45 and randomInt < 0) {
            fail $"RandomInt returned an integer outside the allowed range.";
        }
    }
    
    
    /// # Summary
    /// Checks that @"microsoft.quantum.canon.randomintpow2" obeys ranges.
    operation RandomIntPow2RangeTest () : Unit {
        
        let randIntPow2 = RandomIntPow2(7);
        
        if (randIntPow2 > 127 and randIntPow2 < 0) {
            fail $"RandomIntPow2 returned an integer outside the allowed range.";
        }
    }

    /// # Summary
    /// Checks that @"Microsoft.Quantum.Canon.RandomReal" obeys ranges.
	operation RandomRealTest() : Unit
	{
		for(i in 1..6)
		{
			for(j in 0..10000)
			{
				let numberOfBits = i * 10;
				let random = RandomReal(numberOfBits);
				if(random < 0.0 and random >= 1.0)
				{
					fail $"RandomReal failed with {numberOfBits} bits. Value was {random} which is out of bounds.";
				}
			}
		}
	}
    
}


