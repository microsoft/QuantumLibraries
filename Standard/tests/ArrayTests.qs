// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.
namespace Microsoft.Quantum.Tests {
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Diagnostics;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Logical;
    open Microsoft.Quantum.Math;

    @Test("QuantumSimulator")
    function TestZip() : Unit {

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

    @Test("QuantumSimulator")
    function TestUnzipped() : Unit {
        let first = [6, 5, 5, 3, 2, 1];
        let second = [true, false, false, false, true, false];

        let (first2, second2) = Unzipped(Zip(first, second));
        AllEqualityFactI(first2, first, "Unexpected array of integers");
        AllEqualityFactB(second2, second, "Unexpected array of Booleans");
    }


    @Test("QuantumSimulator")
    function TestLookup() : Unit {

        let array = [1, 12, 71, 103];
        let fn = LookupFunction(array);
        EqualityFactI(fn(0), 1, $"fn(0) did not return array[0]");

        // Make sure we can call in random order!
        EqualityFactI(fn(3), 103, $"fn(3) did not return array[3]");
        EqualityFactI(fn(2), 71, $"fn(2) did not return array[2]");
        EqualityFactI(fn(1), 12, $"fn(1) did not return array[1]");
    }

    internal function AllEqualI(expected : Int[], actual : Int[]) : Bool {
        return All(EqualI, Zip(expected, actual));
    }

    @Test("QuantumSimulator")
    function TestChunks() : Unit {
        let data = [10, 11, 12, 13, 14, 15];

        // 2 × 3 case.
        Fact(All(AllEqualI, Zip(
            [[10, 11], [12, 13], [14, 15]],
            Chunks(2, data)
        )), "Wrong chunks in 2x3 case.");

        // Case with some leftovers.
        Fact(All(AllEqualI, Zip(
            [[10, 11, 12, 13], [14, 15]],
            Chunks(4, data)
        )), "Wrong chunks in case with leftover elements.");
    }

    internal function Squared(x : Int) : Int {
        return x * x;
    }


    function ConstantArrayTest () : Unit {

        let dblArray = ConstantArray(71, 2.17);
        EqualityFactI(Length(dblArray), 71, $"ConstantArray(Int, Double) had the wrong length.");
        let ignore = Mapped(NearEqualityFactD(_, 2.17), dblArray);

        // Stress test by making an array of Int -> Int.
        let fnArray = ConstantArray(7, Squared);
        EqualityFactI(Length(fnArray), 7, $"ConstantArray(Int, Int -> Int) had the wrong length.");
        EqualityFactI(fnArray[3](7), 49, $"ConstantArray(Int, Int -> Int) had the wrong value.");
    }


    function SubarrayTest () : Unit {

        let array0 = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
        let subarrayOdd = Subarray([1, 3, 5, 7, 9], array0);
        let subarrayEven = Subarray([0, 2, 4, 6, 8, 10], array0);
        Fact(All(IsEven, subarrayEven), $"the even elements of [1..10] were not correctly sliced.");
        Fact(not Any(IsEven, subarrayOdd), $"the odd elements of [1..10] were not correctly sliced.");
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

    function SwapOrderToPermuteArrayTest() : Unit {
        let newOrder = [0, 4, 2, 1, 3];
        let expected = [(1, 4), (1, 3)];
        let actual = _SwapOrderToPermuteArray(newOrder);

        EqualityFactI(Length(expected), Length(actual), "Number of swaps does not match");
        for ((exp, act) in Zip(expected, actual)) {
            let (leftExp, rightExp) = exp;
            let (leftAct, rightAct) = act;

            EqualityFactI(leftExp, leftAct, "SWAP doesn't match");
            EqualityFactI(rightExp, rightAct, "SWAP doesn't match");
        }
    }

    function SwappedTest() : Unit {
        let example = [2, 4, 6, 8, 10];
        let expected = [2, 8, 6, 4, 10];
        let leftIndex = 1;
        let rightIndex = 3;
        let newArray = Swapped(leftIndex, rightIndex, example);

        EqualityFactI(Length(expected), Length(newArray), "Swapped array is a different size than original");
        for ((exp, act) in Zip(expected, newArray)) {
            EqualityFactI(exp, act, "Elements did not match");
        }
    }

    function TupleArrayAsNestedArrayTest() : Unit {
        let example = [(0, 1), (2, 3), (4, 5), (6, 7)];
        let expected = [[0, 1], [2, 3], [4, 5], [6, 7]];

        let actual = TupleArrayAsNestedArray(example);
        EqualityFactI(Length(expected), Length(actual), "Arrays are of different sizes");
        for ((exp, act) in Zip(expected, actual)) {
            for ((elementExp, elementAct) in Zip(exp, act)) {
                EqualityFactI(elementExp, elementAct, "Elements did not match");
            }
        }
    }


    function EqualATest() : Unit {
        // arrays of integers
        let equalArrays = EqualA(EqualI, [2, 3, 4], [2, 3, 4]);
        Fact(equalArrays, "Equal arrays were not reported as equal");

        // arrays of doubles
        let differentLength = EqualA(EqualD, [2.0, 3.0, 4.0], [2.0, 3.0]);
        Fact(not differentLength, "Arrays of different length were reported as equal");

        // arrays of Results
        let differentElements = EqualA(EqualR, [One, Zero], [One, One]);
        Fact(not differentElements, "Arrays with different elements were reported as equal");
    }

    @Test("QuantumSimulator")
    operation TestInterleaved() : Unit {
        AllEqualityFactI(Interleaved([1, 2, 3], [-1, -2, -3]), [1, -1, 2, -2, 3, -3], "Interleaving failed");
        AllEqualityFactB(Interleaved(ConstantArray(3, false), ConstantArray(2, true)), [false, true, false, true, false], "Interleaving failed");
    }

    @Test("QuantumSimulator")
    operation TestCumulativeFolded() : Unit {
        AllEqualityFactI(CumulativeFolded(PlusI, 0, SequenceI(1, 5)), [1, 3, 6, 10, 15], "CumulativeFolded failed");
    }

    @Test("QuantumSimulator")
    operation TestTransposed() : Unit {
        for ((actual, expected) in Zip(Transposed([[1, 2, 3], [4, 5, 6]]), [[1, 4], [2, 5], [3, 6]])) {
            AllEqualityFactI(actual, expected, "Transposed failed");
        }
    }

    @Test("QuantumSimulator")
    operation TestColumnVectorAt() : Unit {
        let matrix = [[1, 2, 3], [4, 5, 6], [7, 8, 9]];
        AllEqualityFactI(ColumnVectorAt(0, matrix), [1, 4, 7], "ColumnVectorAt failed");
        AllEqualityFactI(ColumnVectorAt(1, matrix), [2, 5, 8], "ColumnVectorAt failed");
        AllEqualityFactI(ColumnVectorAt(2, matrix), [3, 6, 9], "ColumnVectorAt failed");
    }

    @Test("QuantumSimulator")
    operation TestElementAt() : Unit {
        let lucas = [2, 1, 3, 4, 7, 11, 18, 29, 47, 76];
        let prime = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29];
        let fibonacci = [0, 1, 1, 2, 3, 5, 8, 13, 21, 34];
        let catalan = [1, 1, 2, 5, 14, 42, 132, 429, 1430, 4862];
        let famous2 = Mapped(ElementAt<Int>(2, _), [lucas, prime, fibonacci, catalan]);
        AllEqualityFactI(famous2, [3, 5, 1, 2], "ElementAt failed");
    }

    @Test("QuantumSimulator")
    operation TestElementsAt() : Unit {
        let lucas = [2, 1, 3, 4, 7, 11, 18, 29, 47, 76];
        let prime = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29];
        let fibonacci = [0, 1, 1, 2, 3, 5, 8, 13, 21, 34];
        let catalan = [1, 1, 2, 5, 14, 42, 132, 429, 1430, 4862];
        let famousOdd = Mapped(ElementsAt<Int>(0..2..9, _), [lucas, prime, fibonacci, catalan]);
        for ((actual, expected) in Zip(famousOdd, [[2, 3, 7, 18, 47], [2, 5, 11, 17, 23], [0, 1, 3, 8, 21], [1, 2, 14, 132, 1430]])) {
            AllEqualityFactI(actual, expected, "ElementsAt failed");
        }
    }
}


