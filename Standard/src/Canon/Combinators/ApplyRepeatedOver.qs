// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {
    open Microsoft.Quantum.Diagnostics;
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Intrinsic;
    
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // Helpers to repeatedly apply functions over qubit arrays
    ///////////////////////////////////////////////////////////////////////////////////////////////

    /// # Summary
    /// Applies a list of ops and their targets sequentially on a qubit array.
    ///
    /// # Input
    /// ## listOfOps
    /// List of ops, each taking a qubit array, to be applied. They are applied sequentially, lowest index first.
    /// ## targets
    /// Nested arrays describing the targets of the op. Each array should contain a list of ints describing 
    /// the qubits to be used.
    /// ## register
    /// Qubit register to be acted upon.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApplyOpRepeatedlyOver
    operation ApplySeriesOfOps(listOfOps : (Qubit[] => Unit)[], targets : Int[][], register : Qubit[]) : Unit {
        if (Length(listOfOps) != Length(targets)) {
            fail "The number of ops and number of targets do not match!";
        }
        for (index in 0..Length(listOfOps) - 1) {
            if (Length(targets[index]) > Length(register)) {
                fail "There are too many targets!";
            }
            let opToApply = listOfOps[index];
            let qubitTargets = Subarray(targets[index], register);
            opToApply(qubitTargets);
        }
    }

    /// # Summary
    /// Applies a list of ops and their targets sequentially on a qubit array. (Adjoint)
    ///
    /// # Input
    /// ## listOfOps
    /// List of ops, each taking a qubit array, to be applied. They are applied sequentially, lowest index first.
    /// Each must have a Adjoint functor.
    /// ## targets
    /// Nested arrays describing the targets of the op. Each array should contain a list of ints describing 
    /// the qubits to be used.
    /// ## register
    /// Qubit register to be acted upon.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApplyOpRepeatedlyOver
    operation ApplySeriesOfOpsA(listOfOps : (Qubit[] => Unit is Adj)[], targets : Int[][], register : Qubit[]) : Unit is Adj{
        if (Length(listOfOps) != Length(targets)) {
            fail "The number of ops and number of targets do not match!";
        }
        for (index in 0..Length(listOfOps) - 1) {
            if (Length(targets[index]) > Length(register)) {
                fail "There are too many targets!";
            }
            let opToApply = listOfOps[index];
            let qubitTargets = Subarray(targets[index], register);
            opToApply(qubitTargets);
        }
    }

    /// # Summary
    /// Applies a list of ops and their targets sequentially on a qubit array. (Controlled)
    ///
    /// # Input
    /// ## listOfOps
    /// List of ops, each taking a qubit array, to be applied. They are applied sequentially, lowest index first.
    /// Each must have a Controlled functor.
    /// ## targets
    /// Nested arrays describing the targets of the op. Each array should contain a list of ints describing 
    /// the qubits to be used.
    /// ## register
    /// Qubit register to be acted upon.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApplyOpRepeatedlyOver
    operation ApplySeriesOfOpsC(listOfOps : (Qubit[] => Unit is Ctl)[], targets : Int[][], register : Qubit[]) : Unit is Ctl {
        if (Length(listOfOps) != Length(targets)) {
            fail "The number of ops and number of targets do not match!";
        }
        for (index in 0..Length(listOfOps) - 1) {
            if (Length(targets[index]) > Length(register)) {
                fail "There are too many targets!";
            }
            let opToApply = listOfOps[index];
            let qubitTargets = Subarray(targets[index], register);
            opToApply(qubitTargets);
        }
    }

    /// # Summary
    /// Applies a list of ops and their targets sequentially on a qubit array. (Adjoint + Controlled)
    ///
    /// # Input
    /// ## listOfOps
    /// List of ops, each taking a qubit array, to be applied. They are applied sequentially, lowest index first.
    /// Each must have both a Adjoint and Controlled functor.
    /// ## targets
    /// Nested arrays describing the targets of the op. Each array should contain a list of ints describing 
    /// the qubits to be used.
    /// ## register
    /// Qubit register to be acted upon.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApplyOpRepeatedlyOver
    operation ApplySeriesOfOpsCA(listOfOps : (Qubit[] => Unit is Adj + Ctl)[], targets : Int[][], register : Qubit[]) : Unit is Adj + Ctl {
        if (Length(listOfOps) != Length(targets)) {
            fail "The number of ops and number of targets do not match!";
        }
        for (index in 0..Length(listOfOps) - 1) {
            if (Length(targets[index]) > Length(register)) {
                fail "There are too many targets!";
            }
            let opToApply = listOfOps[index];
            let qubitTargets = Subarray(targets[index], register);
            opToApply(qubitTargets);
        }
    }

    /// # Summary
    /// Applies the same op over a qubit register multiple times.
    ///
    /// # Input
    /// ## op
    /// An operation to be applied multiple times on the qubit register
    /// ## targets
    /// Nested arrays describing the targets of the op. Each array should contain a list of ints describing 
    /// the qubits to be used.
    /// ## register
    /// Qubit register to be acted upon.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApplySeriesOfOps
    operation ApplyOpRepeatedlyOver(op : (Qubit[] => Unit), targets : Int[][], register : Qubit[]) : Unit 
    {
        for (index in 0..Length(targets) - 1) 
        {
            if (Length(targets[index]) > Length(register)) 
            {
                fail "Too many targets!";
            }
            let opTargets = Subarray(targets[index], register);
            op(opTargets);
        }
    }

    /// # Summary
    /// Applies the same op over a qubit register multiple times.
    ///
    /// # Input
    /// ## op
    /// An operation to be applied multiple times on the qubit register
    /// ## targets
    /// Nested arrays describing the targets of the op. Each array should contain a list of ints describing 
    /// the qubits to be used.
    /// ## register
    /// Qubit register to be acted upon.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApplySeriesOfOps
    operation ApplyOpRepeatedlyOverA(op : (Qubit[] => Unit is Adj), targets : Int[][], register : Qubit[]) : Unit is Adj
    {
        for (index in 0..Length(targets) - 1) 
        {
            if (Length(targets[index]) > Length(register)) 
            {
                fail "Too many targets!";
            }
            let opTargets = Subarray(targets[index], register);
            op(opTargets);
        }
    }

    /// # Summary
    /// Applies the same op over a qubit register multiple times.
    ///
    /// # Input
    /// ## op
    /// An operation to be applied multiple times on the qubit register
    /// ## targets
    /// Nested arrays describing the targets of the op. Each array should contain a list of ints describing 
    /// the qubits to be used.
    /// ## register
    /// Qubit register to be acted upon.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApplySeriesOfOps
    operation ApplyOpRepeatedlyOverC(op : (Qubit[] => Unit is Ctl), targets : Int[][], register : Qubit[]) : Unit is Ctl
    {
        for (index in 0..Length(targets) - 1) 
        {
            if (Length(targets[index]) > Length(register)) 
            {
                fail "Too many targets!";
            }
            let opTargets = Subarray(targets[index], register);
            op(opTargets);
        }
    }

    /// # Summary
    /// Applies the same op over a qubit register multiple times.
    ///
    /// # Input
    /// ## op
    /// An operation to be applied multiple times on the qubit register
    /// ## targets
    /// Nested arrays describing the targets of the op. Each array should contain a list of ints describing 
    /// the qubits to be used.
    /// ## register
    /// Qubit register to be acted upon.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApplySeriesOfOps
    operation ApplyOpRepeatedlyOverCA(op : (Qubit[] => Unit is Adj+Ctl), targets : Int[][], register : Qubit[]) : Unit is Adj+Ctl
    {
        for (index in 0..Length(targets) - 1) 
        {
            if (Length(targets[index]) > Length(register)) 
            {
                fail "Too many targets!";
            }
            let opTargets = Subarray(targets[index], register);
            op(opTargets);
        }
    }

    // operation PermuteQubits(targets : Int[][], register : Qubit[]) : Unit {
    //     ApplyOpRepeatedlyOver(ApplyToFirstTwoQubits(SWAP, _), targets, register);
    // }

    /// # Summary
    /// Permutes qubits by using the SWAP operator.
    ///
    /// # Input
    /// ## ordering
    /// Describes the new ordering of the qubits, where the qubit at index i will now be at ordering[i].
    /// ## register
    /// Qubit register to be acted upon.
    operation PermuteQubits(ordering : Int[], register : Qubit[]) : Unit {
        // Do something!
        if (Length(ordering) != Length(register)) {
            fail "The new ordering array has an incorrect number of elements";
        }

        let swapOrder = SwapOrderToPermuteArray(ordering);
        let swapTargets = TupleListToNestedArray(swapOrder);

        let SwapOp = ApplyToFirstTwoQubitsCA(SWAP, _);
        ApplyOpRepeatedlyOverCA(SwapOp, swapTargets, register);
    }

}