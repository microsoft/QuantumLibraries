// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {
    function IntAbs (input : Int) : Int{ 
        mutable tmp = 0;
        if (input < 0) { 
            set tmp = -input; 
        } 
        else {
            set tmp = input; 
        }
        return tmp;             
    }
    
    function IntMax (a : Int, b : Int) : Int { 
        mutable tmp = 0;
        if (a < b) {
            set tmp = b; 
        } 
        else {
            set tmp = a;
        }
        return tmp;
    }
    
    //function IntRange (range : Range) : Int[] {
    //    mutable resultArray = new Int[ IntMax( IntAbs(range.start), IntAbs(range.stop) )]; 
     //   mutable numItems = 0;
     //   for (idx in range) { 
     //       set resultArray[numItems] = idx; 
     //       set numItems = numItems + 1;
     //   }
     //   return resultArray[0..(numItems-1)];
    //}
    
}
