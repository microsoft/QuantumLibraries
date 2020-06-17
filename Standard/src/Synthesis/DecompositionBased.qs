// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.
namespace Microsoft.Quantum.Synthesis {
    open Microsoft.Quantum.Arithmetic;
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Convert;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Logical;
    open Microsoft.Quantum.Math;

    ////////////////////////////////////////////////////////////
    // Decomposition-based synthesis                          //
    ////////////////////////////////////////////////////////////

    internal function DecomposedOn(perm : Int[], var : Int) : ((Int[], Int[]), Int[]) {
        let n = Length(perm);

        mutable visited = new Bool[n];
        mutable left = new Int[n];
        mutable right = new Int[n];

        mutable row = 0;

        while (row != -1) {
            //Message($"Row = {row}");
            //Message($"  set left[{row}] to {(row &&& ~~~(1 <<< var))}");
            set left w/= row <- (row &&& ~~~(1 <<< var));
            set visited w/= row <- true;
            //Message($"  set left[{row ^^^ (1 <<< var)}] to {left[row] ^^^ (1 <<< var)}");
            set left w/= row ^^^ (1 <<< var) <- left[row] ^^^ (1 <<< var);
            set row ^^^= (1 <<< var);
            set visited w/= row <- true;

            //Message($"  set right[{perm[row] ||| (1 <<< var)}] to {perm[row]}");
            set right w/= perm[row] ||| (1 <<< var) <- perm[row];
            //Message($"  set right[{perm[row] &&& ~~~(1 <<< var)}] to {perm[row] ^^^ (1 <<< var)}");
            set right w/= perm[row] &&& ~~~(1 <<< var) <- perm[row] ^^^ (1 <<< var);

            //Message($"  find {perm[row] ^^^ (1 <<< var)} in perm = {perm}   (row = {row})");
            set row = IndexOf(EqualI(perm[row] ^^^ (1 <<< var), _), perm);
            //Message($"  new row = {row}");
            if (visited[row]) {
                set row = IndexOf(EqualB(false, _), visited);
            }
        }

        mutable remainder = new Int[n];
        for ((i, p) in Enumerated(perm)) {
            set remainder w/= left[i] <- right[p];
        }

        return ((left, right), remainder);
    }

    internal function WithZeroInsertedAt (position : Int, vars : Int, x : Int) : Int {
        return ((x &&& (2^(vars-1) - 2^position)) <<< 1) + (x &&& (2^position - 1));
    }

    ////////////////////////////////////////////////////////////
    // Public operation                                       //
    ////////////////////////////////////////////////////////////

    operation ApplyPermutationUsingDecomposition(perm : Int[], qubits : LittleEndian) : Unit {
        let n = BitSizeI(Length(perm) - 1);

        mutable lFunctions = new (BigInt, Int)[n];
        mutable rFunctions = new (BigInt, Int)[n];

        mutable permCopy = perm;
        for (i in 0..n - 1) {
            let ((l, r), remainder) = DecomposedOn(permCopy, i);
            set permCopy = remainder;
            let indices = Mapped(WithZeroInsertedAt(i, n, _), RangeAsIntArray(0..2^(n - 1) - 1));
            set lFunctions w/= i <- (BoolArrayAsBigInt(Mapped(NotEqualI, Subarray(indices, Enumerated(l)))), i);
            set rFunctions w/= i <- (BoolArrayAsBigInt(Mapped(NotEqualI, Subarray(indices, Enumerated(r)))), i);
        }

        let register = qubits!;
        for ((func, target) in lFunctions + Reversed(rFunctions)) {
            if (func != 0L) {
                ControlledXOnTruthTable(func, Exclude([target], register), register[target]);
            }
        }
    }
}
