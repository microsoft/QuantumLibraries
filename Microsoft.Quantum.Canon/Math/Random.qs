// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {
	open Microsoft.Quantum.Primitive;
    open Microsoft.Quantum.Extensions.Convert;
    open Microsoft.Quantum.Extensions.Math;

	// DESIGN NOTES:
	//     Q♭ functions *CANNOT* have side effects. Since random sampling is a side
	//     effect, this immediately implies that any random number generation is modeled
	//     as an operation. That in turn implies that random number generation calls inside
	//     another operation will break adjointability and controllability.

	/// summary:
	///     Generates a random number uniformly sampled in [0, 2^maxBits-1].
	operation RandomIntPow2 (maxBits : Int) : Int
	{
		body {
			mutable number = 0;
			for( idxBit in 0..maxBits - 1 ) {
				let bit = Random([0.5; 0.5]);
				set number = number + bit * 2 ^ idxBit;
			}
			return number;
		}
	}

	/// summary:
	///     Generates a random number uniformly sampled in [0, maxInt).
	operation RandomInt (maxInt : Int) : Int
	{
		body {
			mutable nBits = 0;
			mutable output = 0;
			set nBits = Ceiling(Lg(ToDouble(maxInt)));

			repeat {
				set output = RandomIntPow2(nBits);
			} until(output < maxInt) fixup {}

			return output;
		}

	}
	// Generates a random real number in [0, 1].
	operation RandomReal (bitsRandom : Int) : Double {
		body {
			if (bitsRandom < 1) {
				fail "Number of random bits must be greater than 0.";
			}
			return ToDouble(RandomIntPow2(bitsRandom)) / ToDouble(2 ^ (bitsRandom - 1));
		}
	}

}
