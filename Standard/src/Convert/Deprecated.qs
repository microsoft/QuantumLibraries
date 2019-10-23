// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {
    open Microsoft.Quantum.Convert;

    /// # Deprecated
    /// This function has been removed.
    function AsQubitArray (arr : Qubit[]) : Qubit[] {
        // TODO: add call to _Removed.
        return arr;
    }

    /// # Deprecated
    /// Please use @"microsoft.quantum.convert.resultarrayasint".
    @Deprecated("Microsoft.Quantum.Convert.ResultArrayAsInt")
    function ResultAsInt (results : Result[]) : Int {
        return ResultArrayAsInt(results);
    }

    /// # Deprecated
    /// Please use @"microsoft.quantum.convert.resultasbool".
    @Deprecated("Microsoft.Quantum.Convert.ResultAsBool")
    function BoolFromResult (input : Result) : Bool {
        return ResultAsBool(input);
    }

    /// # Deprecated
    /// Please use @"microsoft.quantum.convert.boolasresult".
    @Deprecated("Microsoft.Quantum.Convert.BoolAsResult")
    function ResultFromBool (input : Bool) : Result {
        return BoolAsResult(input);
    }

    /// # Deprecated
    /// Please use @"microsoft.quantum.convert.intasboolarray".
    @Deprecated("Microsoft.Quantum.Convert.IntAsBoolArray")
    function BoolArrFromPositiveInt (number : Int, bits : Int) : Bool[] {
        return IntAsBoolArray(number, bits);
    }

    /// # Deprecated
    /// Please use @"microsoft.quantum.convert.resultarrayasboolarray".
    @Deprecated("Microsoft.Quantum.Convert.ResultArrayAsBoolArray")
    function BoolArrFromResultArr(input : Result[]) : Bool[] {
        return ResultArrayAsBoolArray(input);
    }

    /// # Deprecated
    /// Please use @"microsoft.quantum.convert.resultarrayasboolarray".
    @Deprecated("Microsoft.Quantum.Convert.BoolArrayAsResultArray")
    function ResultArrFromBoolArr(input : Bool[]) : Result[] {
        return BoolArrayAsResultArray(input);
    }

    /// # Deprecated
    /// Please use @"microsoft.quantum.convert.boolarrayasint".
    @Deprecated("Microsoft.Quantum.Convert.BoolArrayAsInt")
    function PositiveIntFromBoolArr(bits : Bool[]) : Int {
        return BoolArrayAsInt(bits);
    }

    /// # Deprecated
    /// Please use @"microsoft.quantum.canon.resultarrayasint".
    @Deprecated("Microsoft.Quantum.Convert.ResultArrayAsInt")
    function PositiveIntFromResultArr(results : Result[]) : Int {
        return ResultArrayAsInt(results);
    }

    /// # Deprecated
    /// Please use @"microsoft.quantum.convert.boolarrayaspauli".
    @Deprecated("Microsoft.Quantum.Convert.BoolArrayAsPauli")
    function PauliFromBitString(pauli : Pauli, bitApply : Bool, bits : Bool[]) : Pauli[] {
        return BoolArrayAsPauli(pauli, bitApply, bits);
    }

    /// # Deprecated
    /// Please use @"microsoft.quantum.convert.rangeasintarray".
    @Deprecated("Microsoft.Quantum.Convert.RangeAsIntArray")
    function IntArrayFromRange (range: Range) : Int[] {
        return RangeAsIntArray(range);
    }

    /// # Deprecated
    /// Please use @"microsoft.quantum.convert.functionasoperation".
    @Deprecated("Microsoft.Quantum.Convert.FunctionAsOperation")
    function ToOperation<'Input, 'Output>(fn : ('Input -> 'Output)) : ('Input => 'Output) {
        return FunctionAsOperation(fn);
    }

}
