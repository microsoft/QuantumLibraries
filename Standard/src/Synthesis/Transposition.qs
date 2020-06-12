// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.
namespace Microsoft.Quantum.Synthesis {
    open Microsoft.Quantum.Arithmetic;
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Convert;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Logical;


    operation ApplyTransposition(a : Int, b : Int, register : LittleEndian) : Unit is Adj+Ctl {
        if (a != b) {
            let qs = register!;
            let n = Length(qs);

            let abits = IntAsBoolArray(a, n);
            let bbits = IntAsBoolArray(b, n);
            let diff = Where(EqualB(_, true), IntAsBoolArray(a ^^^ b, n));

            within {
                for (targetIndex in Most(diff)) {
                    (ControlledOnBitString(bbits[0..targetIndex - 1] + abits[targetIndex + 1..n - 1], X))(qs[0..targetIndex - 1] + qs[targetIndex + 1..n - 1], qs[targetIndex]);
                }
            } apply {
                let targetIndex = Tail(diff);
                (ControlledOnBitString(bbits[0..targetIndex - 1] + abits[targetIndex + 1..n - 1], X))(qs[0..targetIndex - 1] + qs[targetIndex + 1..n - 1], qs[targetIndex]);
            }
        }
    }
}
