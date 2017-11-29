// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {

    function BoolFromResult( input: Result) : Bool
    {
        if(input == Zero) {
            return false;
        }
        else {
            return true;
        }
    }

    function ResultFromBool( input: Bool) : Result
    {
        if(input == false) {
            return Zero;
        }
        else {
            return One;
        }
    }

    function BoolArrFromResultArr(input : Result[]) : Bool[]
    {
        let nInput = Length(input);
        mutable output = new Bool[nInput];
        for (idx in 0..nInput - 1) {
            set output[idx] = BoolFromResult(input[idx]);
        }
        return output;
    }

    function ResultArrFromBoolArr(input : Bool[]) : Result[]
    {
        let nInput = Length(input);
        mutable output = new Result[nInput];
        for (idx in 0..nInput - 1) {
            set output[idx] = ResultFromBool(input[idx]);
        }
        return output;
    }

    /// # Summary
    /// Produces binary representation of positive integer in little Endian format.
    ///
    /// # Input
    /// ## number
    /// Positive integer.
    /// ## bits
    /// Bits in binary representation of number.
    /// # Remarks
    /// The input "number" must be at most 2^bits -1.
    function BoolArrFromPositiveInt(number : Int, bits : Int) : Bool[]
    {
        //FailOn(number > 2^bits - 1, "Number of output bits must be at least Log_2(integer-1).")

        mutable outputBits = new Bool[bits];
        mutable tempInt = number;

        for ( idxBit in 0..bits - 1 ) {
            if ( tempInt % 2 == 0 ){
                set outputBits[idxBit] = false;
            }
            else {
                set outputBits[idxBit] = true;
            }
            set tempInt = tempInt / 2;
        }

        return outputBits;

    }

    /// # Summary
    /// Produces a positive integer from a string of bits in in little Endian format.
    ///
    /// # Input
    /// ## bits
    /// Bits in binary representation of number.
    function PositiveIntFromBoolArr(bits : Bool[]) : Int
    {
        mutable number = 0;
        let nBits = Length(bits) ;

        for ( idxBit in 0..nBits - 1 ) {
            if (bits[idxBit]) {
                set number = number + 2 ^ idxBit;
            }
        
        }

        return number;

    }

    /// # Summary
    /// Produces a positive integer from a string of Results in in little Endian format.
    ///
    /// # Input
    /// ## results
    /// Results in binary representation of number.
    function PositiveIntFromResultArr(results :Result[]) : Int
    {
        return PositiveIntFromBoolArr(BoolArrFromResultArr(results));
    }
}
