// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.
namespace Microsoft.Quantum.Tests {
    open Microsoft.Quantum.Math;
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Diagnostics;
    
    
    function IsEven (element : Int) : Bool {
        
        return element % 2 == 0;
    }
    
    
    function IsSingleDigit (element : Int) : Bool {
        
        return element >= 0 && element < 10;
    }
    
    
    function Add (input : (Int, Int)) : Int {
        
        let (first, second) = input;
        return first + second;
    }
    
    
    function Squarer (a : Int) : Int {
        
        return a * a;
    }
    
    
    function ForAllTest () : Unit {
        
        ClaimEqualB(All(IsSingleDigit, [3, 4, 7, 8]), true, $"the elements [3, 4, 7, 8] were not found to be single digit numbers.");
        ClaimEqualB(All(IsSingleDigit, [3, 4, 7, 18]), false, $"the elements [3, 4, 7, 18] were found to be single digit numbers.");
    }
    
    
    function ForAnyTest () : Unit {
        
        ClaimEqualB(Any(IsEven, [3, 7, 99, -4]), true, $"the elements [3, 7, 99, -4] were not found to contain at least one even number.");
        ClaimEqualB(Any(IsEven, [3, 7, 99, -41]), false, $"the elements [3, 7, 99, -4] were not found to contain at least one even number.");
    }
    
    
    function FoldTest () : Unit {
        
        let array = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
        ClaimEqualI(Fold(Add, 0, array), 55, $"folding the summation over [1..10] did not yield 55.");
    }
    
    
    function MapTest () : Unit {
        
        let array = [1, 2, 3, 4];
        let squaredArray = Mapped(Squarer, array);
        ClaimEqualI(Fold(Add, 0, squaredArray), 30, $"the sum of the squares of [1, 2, 3, 4] was not found to be 30.");
    }

    function ExtremaTest () : Unit {
        let array = [-10, 10, 7, 0];
        ClaimEqualI(-10, Min(array), $"Min failed.");
        ClaimEqualI(10, Max(array), $"Max failed.");
    }

    function IndexOfTest() : Unit {
        let array = [1, 3, 21, -7, 2, 19];
        let actual = IndexOf(IsEven, array);
        ClaimEqualI(4, actual, $"Expected 4, got {actual}.");
    }

}


