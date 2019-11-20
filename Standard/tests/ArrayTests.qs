// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.
namespace Microsoft.Quantum.Tests {
    open Microsoft.Quantum.Diagnostics;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Arrays;
    
    
    function ZipTest () : Unit {
        
        let left = [1, 2, 101];
        let right = [PauliY, PauliI];
        let zipped = Zip(left, right);
        let (leftActual1, rightActual1) = zipped[0];
        
        if (leftActual1 != 1 or rightActual1 != PauliY) {
            fail $"Expected (1, PauliY), got ({leftActual1}, {rightActual1}).";
        }
        
        let (leftActual2, rightActual2) = zipped[1];
        
        if (leftActual2 != 2 or rightActual2 != PauliI) {
            fail $"Expected (2, PauliI), got ({leftActual2}, {rightActual2}).";
        }
    }
    
    
    function LookupTest () : Unit {
        
        let array = [1, 12, 71, 103];
        let fn = LookupFunction(array);
        EqualityFactI(fn(0), 1, $"fn(0) did not return array[0]");
        
        // Make sure we can call in random order!
        EqualityFactI(fn(3), 103, $"fn(3) did not return array[3]");
        EqualityFactI(fn(2), 71, $"fn(2) did not return array[2]");
        EqualityFactI(fn(1), 12, $"fn(1) did not return array[1]");
    }
    
    
    function ConstantArrayTestHelper (x : Int) : Int {
        
        return x * x;
    }
    
    
    function ConstantArrayTest () : Unit {
        
        let dblArray = ConstantArray(71, 2.17);
        EqualityFactI(Length(dblArray), 71, $"ConstantArray(Int, Double) had the wrong length.");
        let ignore = Mapped(NearEqualityFactD(_, 2.17), dblArray);
        
        // Stress test by making an array of Int -> Int.
        let fnArray = ConstantArray(7, ConstantArrayTestHelper);
        EqualityFactI(Length(fnArray), 7, $"ConstantArray(Int, Int -> Int) had the wrong length.");
        EqualityFactI(fnArray[3](7), 49, $"ConstantArray(Int, Int -> Int) had the wrong value.");
    }
    
    
    function SubarrayTest () : Unit {
        
        let array0 = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
        let subarrayOdd = Subarray([1, 3, 5, 7, 9], array0);
        let subarrayEven = Subarray([0, 2, 4, 6, 8, 10], array0);
        EqualityFactB(All(IsEven, subarrayEven), true, $"the even elements of [1..10] were not correctly sliced.");
        EqualityFactB(Any(IsEven, subarrayOdd), false, $"the odd elements of [1..10] were not correctly sliced.");
        let array1 = [10, 11, 12, 13];
        Ignore(Mapped(EqualityFactI(_, _, $"Subarray failed: subpermutation case."), Zip([12, 11], Subarray([2, 1], array1))));
    }
    
    
    function FilterTest () : Unit {
        
        let array = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
        let evenArray = Filtered(IsEven, array);
        EqualityFactB(All(IsEven, evenArray), true, $"the even elements of [1..10] were not correctly filtered.");
    }
    
    
    function ReverseTest () : Unit {
        
        let array = [1, 2, 3];
        Ignore(Mapped(EqualityFactI(_, _, $"Reversed failed."), Zip([3, 2, 1], Reversed(array))));
    }
    
    
    function ExcludeTest () : Unit {
        
        let array = [10, 11, 12, 13, 14, 15];
        Ignore(Mapped(EqualityFactI(_, _, $"Exclude failed."), Zip([10, 11, 13, 14], Exclude([2, 5], array))));
    }
    
    
    function PadTest () : Unit {
        
        mutable arrayTestCase = [(-5, 2, [10, 11, 12], [10, 11, 12, 2, 2]), (5, 2, [10, 11, 12], [2, 2, 10, 11, 12]), (-3, -2, [10, 11, 12], [10, 11, 12])];
        
        for (idxTest in IndexRange(arrayTestCase)) {
            let (nElementsTotal, defaultElement, inputArray, outputArray) = arrayTestCase[idxTest];
            let paddedArray = Padded(nElementsTotal, defaultElement, inputArray);
            Ignore(Mapped(EqualityFactI(_, _, $"Padded failed."), Zip(outputArray, paddedArray)));
        }
    }

    function EnumeratedTest() : Unit  {
        let example = [37, 12];
        let expected = [(0, 37), (1, 12)];
        let actual = Enumerated(example);

        for ((actualElement, expectedElement) in Zip(actual, expected)) {
            EqualityFactI(Fst(actualElement), Fst(expectedElement), "Indices did not match.");
            EqualityFactI(Snd(actualElement), Snd(expectedElement), "Elements did not match.");
        }
    }

    function SequenceITest() : Unit {
        let example = [(0, 3), (23, 29), (-5, -2)];
        let expected = [[0, 1, 2, 3], [23, 24, 25, 26, 27, 28, 29], [-5, -4, -3, -2]];
        let actual = Mapped(SequenceI, example);

        for ((exp, act) in Zip(expected, actual)) {
            EqualityFactI(Length(exp), Length(act), "Lengths of arrays did not match.");
            for ((i, j) in Zip(exp, act)) {
                EqualityFactI(i, j, "Elements did not match.");
            }
        }
    }

    function SequenceLTest() : Unit {
        let example = [(0L, 3L), (23L, 29L), (-5L, -2L)];
        let expected = [[0L, 1L, 2L, 3L], [23L, 24L, 25L, 26L, 27L, 28L, 29L], [-5L, -4L, -3L, -2L]];
        let actual = Mapped(SequenceL, example);

        for ((exp, act) in Zip(expected, actual)) {
            EqualityFactI(Length(exp), Length(act), "Lengths of arrays did not match.");
            for ((i, j) in Zip(exp, act)) {
                EqualityFactL(i, j, "Elements did not match.");
            }
        }
    }

    function SequenceForNumbersTest() : Unit {
        let example = [3, 5, 0];
        let expected = [[0, 1, 2, 3], [0, 1, 2, 3, 4, 5], [0]];
        let actual = Mapped(SequenceI(0, _), example);

        for ((exp, act) in Zip(expected, actual)) {
            EqualityFactI(Length(exp), Length(act), "Lengths of arrays did not match.");
            for ((i, j) in Zip(exp, act)) {
                EqualityFactI(i, j, "Elements did not match.");
            }
        }
    }

    function IsEmptyTest() : Unit {
        Fact(IsEmpty(new Int[0]), "Empty array marked as non-empty.");
        Fact(IsEmpty(new Qubit[0]), "Empty array marked as non-empty.");
        Fact(IsEmpty(new (Double, (Int -> String))[0]), "Empty array marked as non-empty.");
        Fact(not IsEmpty([PauliX, PauliZ]), "Non-empty array marked as empty.");
        Fact(not IsEmpty([""]), "Non-empty array marked as empty.");
    }

}


