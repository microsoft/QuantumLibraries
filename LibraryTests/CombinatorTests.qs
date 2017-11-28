// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Tests {
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Primitive;

    operation WithTest() : () {
        body {
            let actual = With1(H, X, _);
            let expected = Z;

            AssertOperationsEqualReferenced(ApplyToEach(actual, _), ApplyToEachA(expected, _), 4);
    	}
    }

    // Make sure that if CurryTest fails, it's because of Curry and not
    // something else.
    operation CurryPreTest() : () {
        body {
            AssertOperationsEqualInPlace(Exp([PauliZ], 1.7, _), Exp([PauliZ], 1.7, _), 1);
            AssertOperationsEqualReferenced(Exp([PauliZ], 1.7, _), Exp([PauliZ], 1.7, _), 1);
        }
    }

    operation CurryTest() : () {
        body {
            let curried = CurryOp(Exp([PauliZ], _, _));
            AssertOperationsEqualInPlace(curried(1.7), Exp([PauliZ], 1.7, _), 1);
            AssertOperationsEqualReferenced(curried(1.7), Exp([PauliZ], 1.7, _), 1);
        }
    }

    operation BindTest() : () {
        body {
            let bound = Bind([H; X; H]);
            AssertOperationsEqualReferenced(ApplyToEach(bound, _), ApplyToEachA(Z, _), 3); // FIXME: switched the second argument from ApplyToEach(..) to ApplyToEachA(..)
        }
    }

    operation OperationPowTest() : () {
        body {
            AssertOperationsEqualReferenced(ApplyToEach(OperationPow(H, 2), _), NoOp, 3);
            AssertOperationsEqualReferenced(ApplyToEach(OperationPow(Z, 2), _), NoOp, 3);
            AssertOperationsEqualReferenced(ApplyToEach(OperationPow(S, 4), _), NoOp, 3);
            AssertOperationsEqualReferenced(ApplyToEach(OperationPow(T, 8), _), NoOp, 3);
        }
    }

	function IsEven(element : Int) : Bool {		
		return (element % 2) == 0;
	} 

	function IsSingleDigit(element : Int) : Bool {
		return (element >= 0) && (element < 10);
	}

	function Add(input : (Int, Int)) : Int {
		let (first, second) = input;
		return first + second;
	}

	function Squarer(a: Int) : Int {
		return a * a;
	}

	operation ForAllTest() : () {
		body { 
			AssertBoolEqual(ForAll(IsSingleDigit, [3; 4; 7; 8]), true, "the elements [3; 4; 7; 8] were not found to be single digit numbers.");
			AssertBoolEqual(ForAll(IsSingleDigit, [3; 4; 7; 18]), false, "the elements [3; 4; 7; 18] were found to be single digit numbers.");
		}
	}

	operation ForAnyTest() : () {
		body { 
			AssertBoolEqual(ForAny(IsEven, [3; 7; 99; -4]), true, "the elements [3; 7; 99; -4] were not found to contain at least one even number.");			
			AssertBoolEqual(ForAny(IsEven, [3; 7; 99; -41]), false, "the elements [3; 7; 99; -4] were not found to contain at least one even number.");			
		}
	}

	operation FoldTest() : () {
		body { 
			let array = [1; 2; 3; 4; 5; 6; 7; 8; 9; 10];
		    AssertIntEqual(Fold(Add, 0, array), 55, "folding the summation over [1..10] did not yield 55.");	
		}
	}

	operation SubarrayTest() : () {
		body {
			let array = [0; 1; 2; 3; 4; 5; 6; 7; 8; 9; 10];
			let subarrayOdd = Subarray([1; 3; 5; 7; 9], array); 
			let subarrayEven = Subarray([0; 2; 4; 6; 8; 10], array); 			
			AssertBoolEqual(ForAll(IsEven, subarrayEven), true, "the even elements of [1..10] were not correctly sliced.");
			AssertBoolEqual(ForAny(IsEven, subarrayOdd), false, "the odd elements of [1..10] were not correctly sliced.");
		}
	}

	operation FilterTest() : () {
		body {
			let array = [1; 2; 3; 4; 5; 6; 7; 8; 9; 10];
			let evenArray = Filter(IsEven, array); 
			AssertBoolEqual(ForAll(IsEven, evenArray), true, "the even elements of [1..10] were not correctly filtered.");
		}
	}
	operation MapTest() : () {
		body {
			let array = [1; 2; 3; 4];
			let squaredArray = Map(Squarer, array); 
			AssertIntEqual(Fold(Add, 0, squaredArray), 30, "the sum of the squares of [1; 2; 3; 4] was not found to be 30.");
		}
	}
}
