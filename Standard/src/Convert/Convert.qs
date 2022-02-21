// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Convert {
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Diagnostics;
    open Microsoft.Quantum.Arrays;

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
        return Mapped(ResultAsBool, input);
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
        return Mapped(BoolAsResult, input);
    }

    /// # Summary
    /// Produces a binary representation of a non-negative integer, using the
    /// little-endian representation for the returned array.
    ///
    /// # Input
    /// ## number
    /// A non-negative integer to be converted to an array of boolean values.
    /// ## bits
    /// The number of bits in the binary representation of `number`.
    ///
    /// # Output
    /// An array of boolean values representing `number`.
    ///
    /// # Remarks
    /// The input `bits` must be between 0 and 63.
    /// The input `number` must be between 0 and $2^{\texttt{bits}} - 1$.
    function IntAsBoolArray(number : Int, bits : Int) : Bool[] {
        Fact(bits >= 0 and bits <= 63, $"`bits` must be between 0 and 63 {2^bits}");
        let max = bits < 63 ? 1 <<< bits | 0x7FFFFFFFFFFFFFFF;
        Fact(number >= 0 and number <= max, $"`number` must be between 0 and 2^{bits} - 1, but was {number}.");
        mutable outputBits = [false, size = bits];
        mutable tempInt = number;

        for idxBit in 0 .. bits - 1 {
            set outputBits w/= idxBit <- tempInt % 2 == 0 ? false | true;
            set tempInt = tempInt / 2;
        }

        return outputBits;
    }

    /// # Summary
    /// Produces a non-negative integer from a string of bits in little endian format.
    ///
    /// # Input
    /// ## bits
    /// Bits in binary representation of number.
    function BoolArrayAsInt(bits : Bool[]) : Int {
        Fact(Length(bits) < 64, $"`Length(bits)` must be less than 64, but was {Length(bits)}.");
        mutable number = 0;
        let nBits = Length(bits);

        for idxBit in 0 .. nBits - 1 {
            if (bits[idxBit]) {
                set number = number + 2 ^ idxBit;
            }
        }

        return number;
    }

    /// # Summary
    /// Produces a non-negative integer from a string of Results in little endian format.
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
    function FunctionAsOperation<'Input, 'Output>(fn : ('Input -> 'Output)) : ('Input => 'Output) {
        return Call(fn, _);
    }

    /// # Summary
    /// Given a bit string, returns a multi-qubit Pauli operator
    /// represented as an array of single-qubit Pauli operators.
    ///
    /// # Input
    /// ## pauli
    /// Pauli operator to apply to qubits where `bitsApply == bits[idx]`.
    /// ## bitApply
    /// apply Pauli if bit is this value.
    /// ## bits
    /// Boolean array.
    ///
    /// # Remarks
    /// The Boolean array and the quantum register must be of equal length.
    function BoolArrayAsPauli(pauli : Pauli, bitApply : Bool, bits : Bool[]) : Pauli[] {
        let nBits = Length(bits);
        mutable paulis = [PauliI, size = nBits];

        for idxBit in 0 .. nBits - 1 {
            set paulis w/= idxBit <- bits[idxBit] == bitApply ? pauli | PauliI;
        }

        return paulis;
    }

    /// # Summary
    /// Creates an array `arr` of integers enumerated by start..step..end.
    ///
    /// # Input
    /// ## range
    /// A `Range` of values `start..step..end` to be converted to an array.
    ///
    /// # Output
    /// A new array of integers corresponding to values iterated over by `range`.
    ///
    /// # Remarks
    /// ## Example
    /// ```qsharp
    /// // The following returns [1,3,5,7];
    /// let array = RangeAsIntArray(1..2..8);
    /// ```
    function RangeAsIntArray(range : Range) : Int[] {
        let start = RangeStart(range);
        let step = RangeStep(range);
        let end = RangeEnd(range);
        if (IntAsDouble(end - start) / IntAsDouble(step) >= 0.0) {
            let nTerms = (end - start) / step + 1;
            mutable array = [0, size = nTerms];
            for idx in 0..nTerms - 1 {
               set array w/= idx <- start + idx * step;
            }
            return array;
        }
        else {
            return [0, size = 0];
        }
    }

    /// # Summary
    /// Converts a given integer number to an equivalent string representation.
    ///
    /// # Description
    /// Returns a string given a BigInt.
    ///
    /// # Input
    /// ## a
    /// The big integer to be represented as a string.
    ///
    /// # Output
    /// The value of `a` formatted as a string.
    ///
    /// # Example
    /// ```qsharp
    /// let nAsString = BigIntAsString(12345678901234567890L);
    /// // Displays 12345678901234567890.
    /// Message(nAsString);
    /// ```
    function BigIntAsString(a : BigInt) : String {
         return $"{a}";
    }

}
