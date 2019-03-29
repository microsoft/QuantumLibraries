// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {

    /// # Deprecated
    /// This function has been removed.
    function AsQubitArray (arr : Qubit[]) : Qubit[] {
        // TODO: add call to _Removed.
        return arr;
    }

    /// # Deprecated
    /// Please use @"microsoft.quantum.canon.resultarrayasint".
    function ResultAsInt (results : Result[]) : Int {
        Renamed("Microsoft.Quantum.Canon.ResultAsInt", "Microsoft.Quantum.Canon.ResultArrayAsInt");
        return ResultArrayAsInt(results);
    }

    /// # Deprecated
    /// Please use @"microsoft.quantum.canon.resultasbool".
    function BoolFromResult (input : Result) : Bool {
        Renamed("Microsoft.Quantum.Canon.BoolFromResult", "Microsoft.Quantum.Canon.ResultAsBool");
        return ResultAsBool(input);
    }

    /// # Deprecated
    /// Please use @"microsoft.quantum.canon.boolasresult".
    function ResultFromBool (input : Bool) : Result {
        Renamed("Microsoft.Quantum.Canon.ResultFromBool", "Microsoft.Quantum.Canon.BoolAsResult");
        return BoolAsResult(input);
    }

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Canon.IntAsBoolArray".
    function BoolArrFromPositiveInt (number : Int, bits : Int) : Bool[] {
        Renamed("Microsoft.Quantum.Canon.BoolArrFromPositiveInt", "Microsoft.Quantum.Canon.IntAsBoolArray");
        return IntAsBoolArray(number, bits);
    }

}
