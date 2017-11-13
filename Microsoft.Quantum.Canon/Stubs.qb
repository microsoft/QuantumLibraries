// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {
	open Microsoft.Quantum.Primitive;
    open Microsoft.Quantum.Extensions.Convert;
    open Microsoft.Quantum.Extensions.Math;

	// FIXME: make this a function!
	// Perhaps ToDouble as Float usually refers to single-precision, and we cannot overload Double()?
	function Float(value : Int)  : Double
	{
		return ToDouble(value);
	}

	operation MultiM(targets : Qubit[]) : Result[]
	{
		body{
			mutable results = new Result[Length(targets)];
			for(idxQubit in 0..Length(targets)-1){
				set results[idxQubit] = M(targets[idxQubit]);
			}
			return results;
		}
	}

    operation HY(target : Qubit) : () {
        body {
            H(target);
            S(target);
        }

        adjoint auto
        controlled auto
        controlled adjoint auto
    }

}
