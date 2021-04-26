// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.
namespace Microsoft.Quantum.Tests {
    open Microsoft.Quantum.Math;
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Diagnostics;
    
    
    function IsEven (element : Int) : Bool {
        
        return element % 2 == 0;
    }
    
    
    function IsSingleDigit (element : Int) : Bool {
        
        return element >= 0 and element < 10;
    }
    
    
    function Add (input : (Int, Int)) : Int {
        
        let (first, second) = input;
        return first + second;
    }
    
    
    function Squarer (a : Int) : Int {
        
        return a * a;
    }
    
    @Test("QuantumSimulator")
    function ForAllIsCorrect() : Unit {
        
        EqualityFactB(All(IsSingleDigit, [3, 4, 7, 8]), true, $"the elements [3, 4, 7, 8] were not found to be single digit numbers.");
        EqualityFactB(All(IsSingleDigit, [3, 4, 7, 18]), false, $"the elements [3, 4, 7, 18] were found to be single digit numbers.");
    }
        
    @Test("QuantumSimulator")
    function ForAnyIsCorrect() : Unit {
        
        EqualityFactB(Any(IsEven, [3, 7, 99, -4]), true, $"the elements [3, 7, 99, -4] were not found to contain at least one even number.");
        EqualityFactB(Any(IsEven, [3, 7, 99, -41]), false, $"the elements [3, 7, 99, -4] were not found to contain at least one even number.");
    }
    
    @Test("QuantumSimulator")
    function FoldIsCorrect() : Unit {
        
        let array = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
        EqualityFactI(Fold(Add, 0, array), 55, $"folding the summation over [1..10] did not yield 55.");
    }
    
    @Test("QuantumSimulator")
    function MapIsCorrect() : Unit {
        
        let array = [1, 2, 3, 4];
        let squaredArray = Mapped(Squarer, array);
        EqualityFactI(Fold(Add, 0, squaredArray), 30, $"the sum of the squares of [1, 2, 3, 4] was not found to be 30.");
    }

    @Test("QuantumSimulator")
    function TestMappedOverNonEmptyRange() : Unit {
        AllEqualityFactI(MappedOverRange(PlusI(_, 2), 1..5), [3, 4, 5, 6, 7], "MappedOverRange failed.");
    }

    @Test("QuantumSimulator")
    function TestMappedOverReversedRange() : Unit {
        AllEqualityFactI(MappedOverRange(TimesI(_, 2), 4..-2..-4), [8, 4, 0, -4, -8], "MappedOverRange failed.");
    }

    @Test("QuantumSimulator")
    function TestMappedOverEmpty() : Unit {
        AllEqualityFactI(MappedOverRange(TimesI(_, 2), 1..-1..2), new Int[0], "MappedOverRange failed.");
    }

    @Test("QuantumSimulator")
    function TestFlatMapped() : Unit {
        let numbers = FlatMapped(SequenceI(1, _), SequenceI(1, 5));
        AllEqualityFactI(numbers, [1, 1, 2, 1, 2, 3, 1, 2, 3, 4, 1, 2, 3, 4, 5], "FlatMapped failed");
    }

    @Test("QuantumSimulator")
    function TestFlattened() : Unit {
        let numbers = Flattened(ConstantArray(3, SequenceI(1, 3)));
        AllEqualityFactI(numbers, [1, 2, 3, 1, 2, 3, 1, 2, 3], "Flattened failed");
    }

    @Test("QuantumSimulator")
    function ExtremaIsCorrect() : Unit {
        let array = [-10, 10, 7, 0];
        EqualityFactI(-10, Min(array), $"Min failed.");
        EqualityFactI(10, Max(array), $"Max failed.");
    }

    @Test("QuantumSimulator")
    function IndexOfIsCorrect() : Unit {
        let array = [1, 3, 21, -7, 2, 19];
        let actual = IndexOf(IsEven, array);
        EqualityFactI(4, actual, $"Expected 4, got {actual}.");
    }

}


