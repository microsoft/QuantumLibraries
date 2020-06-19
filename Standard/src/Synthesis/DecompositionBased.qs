// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Synthesis {
    open Microsoft.Quantum.Arithmetic;
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Convert;
    open Microsoft.Quantum.Diagnostics;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Logical;
    open Microsoft.Quantum.Math;

    ////////////////////////////////////////////////////////////
    // Decomposition-based synthesis                          //
    ////////////////////////////////////////////////////////////

    /// # Summary
    /// Decomposes a permutation on a variable
    ///
    /// # Description
    /// Given a permutation $\pi$ (`perm`) and an index $i$ (`index`), this method
    /// returns three permutations $((\pi_l, \pi_r), \pi')$ such that the images
    /// of $\pi_l$ and $\pi_r$ do not change bits in their elements at indexes other
    /// than $i$ and images of $\pi'$ do not change bit $i$ in their elements.
    internal function DecomposedOn(perm : Int[], index : Int) : ((Int[], Int[]), Int[]) {
        let n = Length(perm);

        mutable visited = new Bool[n];
        mutable left = new Int[n];
        mutable right = new Int[n];

        mutable row = 0;

        while (row != -1) {
            set left w/= row <- (row &&& ~~~(1 <<< index));
            set visited w/= row <- true;
            set left w/= row ^^^ (1 <<< index) <- left[row] ^^^ (1 <<< index);
            set row ^^^= (1 <<< index);
            set visited w/= row <- true;

            set right w/= perm[row] ||| (1 <<< index) <- perm[row];
            set right w/= perm[row] &&& ~~~(1 <<< index) <- perm[row] ^^^ (1 <<< index);

            set row = IndexOf(EqualI(perm[row] ^^^ (1 <<< index), _), perm);
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

    /// # Summary
    /// State during decomposition based on variable indexes
    ///
    /// # Description
    /// The state holds the current permutation and the currently generated functions
    /// for controlled gates on the left, and controlled gates on the right.
    internal newtype DecompositionState = (perm : Int[], lfunctions : (BigInt, Int)[], rfunctions : (BigInt, Int)[]);

    /// # Summary
    /// Decomposition logic for a single variable index
    ///
    /// # Description
    /// This takes the current state and generates an updated permutation
    /// and possibly adds new functions for controlled gates.
    internal function GetTruthTablesForGatesFolder(numVars : Int, state : DecompositionState, index : Int) : DecompositionState {
        let ((l, r), remainder) = DecomposedOn(state::perm, index);

        let indices = Mapped(WithZeroInsertedAt(index, _), SequenceI(0, 2^(numVars - 1) - 1));
        let lFunc = BoolArrayAsBigInt(Mapped(NotEqualI, Subarray(indices, Enumerated(l))) + [false]); // make sure that lFunc is positive
        let rFunc = BoolArrayAsBigInt(Mapped(NotEqualI, Subarray(indices, Enumerated(r))) + [false]); // make sure that lFunc is positive

        return DecompositionState(
            remainder,
            lFunc == 0L ? state::lfunctions | state::lfunctions + [(lFunc, index)],
            rFunc == 0L ? state::rfunctions | [(rFunc, index)] + state::rfunctions
        );
    }

    /// # Summary
    /// Collect all functions for controlled gates by folding through all variable indexes
    internal function GetTruthTablesForGates (perm : Int[], variableOrder : Int[]) : (BigInt, Int)[] {
        let numVars = Length(variableOrder);

        let initialState = DecompositionState(perm, new (BigInt, Int)[0], new (BigInt, Int)[0]);
        let (_, lfunctions, rfunctions) = (Fold(GetTruthTablesForGatesFolder(numVars, _, _), initialState, variableOrder))!;

        return lfunctions + rfunctions;
    }

    ////////////////////////////////////////////////////////////
    // Public operation                                       //
    ////////////////////////////////////////////////////////////

    /// # Summary
    /// Permutes the amplitudes in a quantum state given a permutation
    /// using decomposition-based synthesis.
    ///
    /// # Description
    /// This procedure implements the decomposition based
    /// synthesis approach.  Input is a permutation $\pi$ over $2^n$ elements
    /// $\{0, \dots, 2^n-1\}$, which represents an $n$-variable reversible Boolean function. 
    /// The algorithm performs iteratively the following steps for each variable
    /// index $i$:
    ///
    /// 1. Compute $((\pi_l, \pi_r), \pi')$ such that the images
    ///    of $\pi_l$ and $\pi_r$ do not change bits in their elements at indexes other
    ///    than $i$ and images of $\pi'$ do not change bit $i$ in their elements.
    /// 2. Set $\pi \leftarrow \pi'$, and derive truth tables from $\pi_l$ and $\pi_r$
    ///    based on elements that are not fix-points.
    ///
    /// After applying these steps for all variable indexes, the remaining
    /// permutation `\pi` will be the identity, and based on the collectd truth
    /// tables and indexes, one can apply truth-table controlled @"microsoft.quantum.intrinsic.x"
    /// operations using the @"microsoft.quantum.synthesis.applyxcontrolledontruthtable" operation.
    ///
    /// The variable order is $0, \dots, n - 1$.  A custom variable order can be specified
    /// in the operation @"microsoft.quantum.synthesis.applypermutationusingdecompositionwithvariableorder".
    ///
    /// # Input
    /// ## perm
    /// A permutation of $2^n$ elements starting from 0.
    /// ## qubits
    /// A list of $n$ qubits to which the permutation is applied to.
    ///
    /// # Example
    /// To synthesize a `SWAP` operation:
    /// ```Q#
    /// using (qubits = Qubit[2]) {
    ///   ApplyPermutationUsingDecomposition([0, 2, 1, 3], qubits);
    /// }
    /// ```
    ///
    /// # References
    /// - [*Alexis De Vos*, *Yvan Van Rentergem*,
    ///    Adv. in Math. of Comm. 2(2), 2008, pp. 183--200](http://www.aimsciences.org/article/doi/10.3934/amc.2008.2.183)
    /// - [*Mathias Soeken*, *Laura Tague*, *Gerhard W. Dueck*, *Rolf Drechsler*,
    ///    Journal of Symbolic Computation 73 (2016), pp. 1--26](https://www.sciencedirect.com/science/article/pii/S0747717115000188?via%3Dihub)
    ///
    /// # See Also
    /// - Microsoft.Quantum.Synthesis.ApplyPermutationUsingDecompositionWithVariableOrder
    operation ApplyPermutationUsingDecomposition(perm : Int[], qubits : LittleEndian) : Unit is Adj+Ctl {
        ApplyPermutationUsingDecompositionWithVariableOrder(perm, SequenceI(0, Length(qubits!) - 1), qubits);
    }

    /// # Summary
    /// Permutes the amplitudes in a quantum state given a permutation
    /// using decomposition-based synthesis.
    ///
    /// # Description
    /// This operation is a more general version of @"microsoft.quantum.synthesis.applypermutationusingdecomposition"
    /// in which the variable order can be specified.
    ///
    /// # Input
    /// ## perm
    /// A permutation of $2^n$ elements starting from 0.
    /// ## variableOrder
    /// A permutation of $n$ elements starting from 0.
    /// ## qubits
    /// A list of $n$ qubits to which the permutation is applied to.
    ///
    /// # Example
    /// To synthesize a `SWAP` operation:
    /// ```Q#
    /// using (qubits = Qubit[2]) {
    ///   ApplyPermutationUsingDecompositionWithVariableOrder([0, 2, 1, 3], [1, 0], qubits);
    /// }
    /// ```
    ///
    /// # See Also
    /// - Microsoft.Quantum.Synthesis.ApplyPermutationUsingDecomposition
    operation ApplyPermutationUsingDecompositionWithVariableOrder(perm : Int[], variableOrder : Int[], qubits : LittleEndian) : Unit is Adj+Ctl {
        Fact(IsPermutation(perm), "perm must be a permutation");
        EqualityFactI(Length(perm), 2^Length(qubits!), $"Length of perm must be {2^Length(qubits!)}");

        Fact(IsPermutation(variableOrder), "variableOrder must be a permutation");
        EqualityFactI(Length(variableOrder), Length(qubits!), $"Length of variableOrder must be {Length(qubits!)}");

        let register = qubits!;

        for ((func, target) in GetTruthTablesForGates(perm, variableOrder)) {
            ApplyXControlledOnTruthTable(func, Exclude([target], register), register[target]);
        }
    }
}
