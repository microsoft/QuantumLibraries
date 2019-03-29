// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {

    /// # Summary
    /// Converts a `Result` type to a `Bool` type, where `One` is mapped to
    /// `true` and `Zero` is mapped to `false`.
    ///
    /// # Input
    /// ## input
    /// `Result` to be converted.
    ///
    /// # Output
    /// A `Bool` representing the `input`.
    function ResultAsBool(input : Result) : Bool {
        return (input == Zero ? false | true);
    }

    /// # Summary
    /// Converts a `Bool` type to a `Result` type, where `true` is mapped to
    /// `One` and `false` is mapped to `Zero`.
    ///
    /// # Input
    /// ## input
    /// `Bool` to be converted.
    ///
    /// # Output
    /// A `Result` representing the `input`.
    function BoolAsResult (input : Bool) : Result {
        return (input == false ? Zero | One);
    }

    /// # Summary
    /// Converts a `Result[]` type to a `Bool[]` type, where `One` is mapped to
    /// `true` and `Zero` is mapped to `false`.
    ///
    /// # Input
    /// ## input
    /// `Result[]` to be converted.
    ///
    /// # Output
    /// A `Bool[]` representing the `input`.
    function ResultArrayAsBoolArray(input : Result[]) : Bool[] {
        return Map(ResultAsBool, input);
    }

    /// # Summary
    /// Converts a `Bool[]` type to a `Result[]` type, where `true` is mapped to
    /// `One` and `false` is mapped to `Zero`.
    ///
    /// # Input
    /// ## input
    /// `Bool[]` to be converted.
    ///
    /// # Output
    /// A `Result[]` representing the `input`.
    function BoolArrayAsResultArray(input : Bool[]) : Result[] {
        return Map(BoolAsResult, input);
    }

    /// # Summary
    /// Produces a binary representation of a positive integer, using the
    /// little-endian representation for the returned array.
    ///
    /// # Input
    /// ## number
    /// A positive integer to be converted to an array of boolean values.
    /// ## bits
    /// The number of bits in the binary representation of `number`.
    ///
    /// # Output
    /// An array of boolean values representing `number`.
    ///
    /// # Remarks
    /// The input `number` must be at most $2^{\texttt{bits}} - 1$.
    function IntAsBoolArray(number : Int, bits : Int) : Bool[] {
        AssertBoolEqual(number >= 0 && number < 2 ^ bits, true, $"`number` must be between 0 and 2^`bits` - 1");
        mutable outputBits = new Bool[bits];
        mutable tempInt = number;

        for (idxBit in 0 .. bits - 1) {
            set outputBits[idxBit] = tempInt % 2 == 0 ? false | true;
            set tempInt = tempInt / 2;
        }

        return outputBits;
    }

    /// # Summary
    /// Produces a positive integer from a string of bits in little endian format.
    ///
    /// # Input
    /// ## bits
    /// Bits in binary representation of number.
    function BoolArrayAsInt(bits : Bool[]) : Int {
        AssertBoolEqual(Length(bits) < 64, true, $"`Length(bits)` must be less than 64");
        mutable number = 0;
        let nBits = Length(bits);

        for (idxBit in 0 .. nBits - 1) {
            if (bits[idxBit]) {
                set number = number + 2 ^ idxBit;
            }
        }

        return number;
    }

    /// # Summary
    /// Produces a positive integer from a string of Results in little endian format.
    ///
    /// # Input
    /// ## results
    /// Results in binary representation of number.
    function ResultArrayAsInt(results : Result[]) : Int {
        return BoolArrayAsInt(ResultArrayAsBoolArray(results));
    }

    /// # Summary
    /// Calls a function with a given input.
    ///
    /// # Description
    /// Given a function and an input to that function, calls the function
    /// and returns its output.
    ///
    /// # Input
    /// ## fn
    /// A function to be called.
    /// ## input
    /// The input to be passed to the function.
    ///
    /// # Output
    /// The result of calling `fn`.
    ///
    /// # Remarks
    /// This operation is mainly useful for forcing a function to be called
    /// at a specific place within an operation, or for calling a function
    /// where an operation is expected.
    operation Call<'Input, 'Output>(fn : ('Input -> 'Output), input : 'Input) : 'Output {
        return fn(input);
    }

    /// # Summary
    /// Converts functions to operations.
    ///
    /// # Description
    /// Given a function, returns an operation which calls that function,
    /// and which does nothing else.
    ///
    /// # Input
    /// ## fn
    /// A function to be converted to an operation.
    ///
    /// # Output
    /// An operation `op` such that `op(input)` is identical to `fn(input)`
    /// for all `input`.
    ///
    /// # Type Parameters
    /// ## 'Input
    /// Input type of the function to be converted.
    /// ## 'Output
    /// Output type of the function to be converted.
    ///
    /// # Remarks
    /// This is mainly useful for passing functions to functions or operations
    /// which expect an operation as input.
    function ToOperation<'Input, 'Output>(fn : ('Input -> 'Output)) : ('Input => 'Output) {
        return Call(fn, _);
    }

}
