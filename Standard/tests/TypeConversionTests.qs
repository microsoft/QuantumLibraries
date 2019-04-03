// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.
namespace Microsoft.Quantum.Tests {
    open Microsoft.Quantum.Primitive;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Convert;
    open Microsoft.Quantum.Diagnostics;

    function Square(x : Int) : Int {
        return x * x;
    }

    operation ApplyOp<'T, 'U>(op : ('T => 'U), input : 'T) : 'U {
        return op(input);
    }

    operation CallTest() : Unit {
        ClaimEqualI(Call(Square, 4), 16, "Call failed with Square.");
    }

    operation ToOperationTest() : Unit {
        let op = FunctionAsOperation(Square);
        ClaimEqualI(ApplyOp(op, 3), 9, "ToOperation failed with Square.");
    }

}


