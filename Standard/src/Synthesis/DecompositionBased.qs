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
            set left w/= row <- (row &&& ~~~(1 <<< var));
            set visited w/= row <- true;
            set left w/= row ^^^ (1 <<< var) <- left[row] ^^^ (1 <<< var);
            set row ^^^= (1 <<< var);
            set visited w/= row <- true;

            set right w/= perm[row] ||| (1 <<< var) <- perm[row];
            set right w/= perm[row] &&& ~~~(1 <<< var) <- perm[row] ^^^ (1 <<< var);

            set row = IndexOf(EqualI(perm[row] ^^^ (1 <<< var), _), perm);
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

    /// # Summary
    /// Insert a 0-bit into an integer
    ///
    /// # Description
    /// This operation takes an integer, inserts a 0 at bit `position`, and returns
    /// the updated value as an integer.  For example, inserting a 0 at position 2
    /// in the number 10 ($10_{10} = 1010_{2}$) returns the number 18 ($18_{10} = 10010_{2}$).
    ///
    /// # Input
    /// ## position
    /// The position at which 0 is inserted
    /// ## value
    /// The value that is modified
    ///
    /// # Output
    /// Modified value
    internal function WithZeroInsertedAt (position : Int, value : Int) : Int {
        let mask = 2^position - 1;
        return ((value &&& ~~~mask) <<< 1) + (value &&& mask);
    }

    internal function GetTruthTablesForGates (perm : Int[]) : (BigInt, Int)[] {
        let n = BitSizeI(Length(perm) - 1);

        mutable lFunctions = new (BigInt, Int)[n];
        mutable rFunctions = new (BigInt, Int)[n];

        mutable permCopy = perm;
        for (i in 0..n - 1) {
            let ((l, r), remainder) = DecomposedOn(permCopy, i);
            set permCopy = remainder;
            let indices = Mapped(WithZeroInsertedAt(i, _), SequenceI(0, 2^(n - 1) - 1));
            
            let lFunc = BoolArrayAsBigInt(Mapped(NotEqualI, Subarray(indices, Enumerated(l))));
            let rFunc = BoolArrayAsBigInt(Mapped(NotEqualI, Subarray(indices, Enumerated(r))));

            if (lFunc != 0L) {
                set lFunctions w/= i <- (lFunc, i);
            }
            if (rFunc != 0L) {
                set rFunctions w/= i <- (rFunc, i);
            }
        }

        return lFunctions + Reversed(rFunctions);
    }

    ////////////////////////////////////////////////////////////
    // Public operation                                       //
    ////////////////////////////////////////////////////////////

    operation ApplyPermutationUsingDecomposition(perm : Int[], qubits : LittleEndian) : Unit is Adj+Ctl {
        let register = qubits!;
        for ((func, target) in GetTruthTablesForGates(perm)) {
            ApplyXControlledOnTruthTable(func, Exclude([target], register), register[target]);
        }
    }
}
