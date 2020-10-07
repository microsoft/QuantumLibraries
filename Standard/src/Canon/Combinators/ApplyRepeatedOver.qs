// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {
    open Microsoft.Quantum.Diagnostics;
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Intrinsic;
    
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // Helpers to repeatedly apply operations over qubit arrays
    ///////////////////////////////////////////////////////////////////////////////////////////////

    /// # Summary
    /// Applies a list of ops and their targets sequentially on an array.
    ///
    /// # Input
    /// ## listOfOps
    /// List of ops, each taking a 'T array, to be applied. They are applied sequentially, lowest index first.
    /// ## targets
    /// Nested arrays describing the targets of the op. Each array should contain a list of ints describing 
    /// the indices to be used.
    /// ## register
    /// Qubit register to be acted upon.
    ///
    /// ## Example
    /// // The following applies Exp([PauliX, PauliY], 0.5) to qubits 0, 1
    /// // then X to qubit 2
    /// let ops = [Exp([PauliX, PauliY], 0.5, _), ApplyToFirstQubit(X, _)];
    /// let indices = [[0, 1], [2]];
    /// ApplySeriesOfOps(ops, indices, qubitArray);
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApplyOpRepeatedlyOver
    operation ApplySeriesOfOps<'T>(listOfOps : ('T[] => Unit)[], targets : Int[][], register : 'T[]) : Unit {
        if (Length(listOfOps) != Length(targets)) {
            fail "The number of ops and number of targets do not match!";
        }
        for ((op, targetIndices) in Zipped(listOfOps, targets)) {
            if (Length(targetIndices) > Length(register)) {
                fail "There are too many targets!";
            }
            op(Subarray(targetIndices, register));
        }
    }

    /// # Summary
    /// Applies a list of ops and their targets sequentially on an array. (Adjoint)
    ///
    /// # Input
    /// ## listOfOps
    /// List of ops, each taking a 'T array, to be applied. They are applied sequentially, lowest index first.
    /// Each must have an adjoint functor
    /// ## targets
    /// Nested arrays describing the targets of the op. Each array should contain a list of ints describing 
    /// the indices to be used.
    /// ## register
    /// Qubit register to be acted upon.
    ///
    /// ## Example
    /// // The following applies Exp([PauliX, PauliY], 0.5) to qubits 0, 1
    /// // then X to qubit 2
    /// let ops = [Exp([PauliX, PauliY], 0.5, _), ApplyToFirstQubitA(X, _)];
    /// let indices = [[0, 1], [2]];
    /// ApplySeriesOfOpsA(ops, indices, qubitArray);
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApplyOpRepeatedlyOver
    operation ApplySeriesOfOpsA<'T>(listOfOps : ('T[] => Unit is Adj)[], targets : Int[][], register : 'T[]) : Unit is Adj{
        if (Length(listOfOps) != Length(targets)) {
            fail "The number of ops and number of targets do not match!";
        }
        for ((op, targetIndices) in Zipped(listOfOps, targets)) {
            if (Length(targetIndices) > Length(register)) {
                fail "There are too many targets!";
            }
            op(Subarray(targetIndices, register));
        }
    }

    /// # Summary
    /// Applies a list of ops and their targets sequentially on an array. (Controlled)
    ///
    /// # Input
    /// ## listOfOps
    /// List of ops, each taking a 'T array, to be applied. They are applied sequentially, lowest index first.
    /// Each must have a Controlled functor
    /// ## targets
    /// Nested arrays describing the targets of the op. Each array should contain a list of ints describing 
    /// the indices to be used.
    /// ## register
    /// Qubit register to be acted upon.
    ///
    /// ## Example
    /// // The following applies Exp([PauliX, PauliY], 0.5) to qubits 0, 1
    /// // then X to qubit 2
    /// let ops = [Exp([PauliX, PauliY], 0.5, _), ApplyToFirstQubitC(X, _)];
    /// let indices = [[0, 1], [2]];
    /// ApplySeriesOfOpsC(ops, indices, qubitArray);
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApplyOpRepeatedlyOver
    operation ApplySeriesOfOpsC<'T>(listOfOps : ('T[] => Unit is Ctl)[], targets : Int[][], register : 'T[]) : Unit is Ctl{
        if (Length(listOfOps) != Length(targets)) {
            fail "The number of ops and number of targets do not match!";
        }
        for ((op, targetIndices) in Zipped(listOfOps, targets)) {
            if (Length(targetIndices) > Length(register)) {
                fail "There are too many targets!";
            }
            op(Subarray(targetIndices, register));
        }
    }

    /// # Summary
    /// Applies a list of ops and their targets sequentially on an array. (Adjoint + Controlled)
    ///
    /// # Input
    /// ## listOfOps
    /// List of ops, each taking a 'T array, to be applied. They are applied sequentially, lowest index first.
    /// Each must have both an Adjoint and Controlled functor.
    /// ## targets
    /// Nested arrays describing the targets of the op. Each array should contain a list of ints describing 
    /// the indices to be used.
    /// ## register
    /// Qubit register to be acted upon.
    ///
    /// ## Example
    /// // The following applies Exp([PauliX, PauliY], 0.5) to qubits 0, 1
    /// // then X to qubit 2
    /// let ops = [Exp([PauliX, PauliY], 0.5, _), ApplyToFirstQubitCA(X, _)];
    /// let indices = [[0, 1], [2]];
    /// ApplySeriesOfOpsCA(ops, indices, qubitArray);
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApplyOpRepeatedlyOver
    operation ApplySeriesOfOpsCA<'T>(listOfOps : ('T[] => Unit is Adj + Ctl)[], targets : Int[][], register : 'T[]) : Unit is Adj + Ctl{
        if (Length(listOfOps) != Length(targets)) {
            fail "The number of ops and number of targets do not match!";
        }
        for ((op, targetIndices) in Zipped(listOfOps, targets)) {
            if (Length(targetIndices) > Length(register)) {
                fail "There are too many targets!";
            }
            op(Subarray(targetIndices, register));
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
        for (target in targets) 
        {
            if (Length(target) > Length(register)) 
            {
                fail "Too many targets!";
            }
            let opTargets = Subarray(target, register);
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
        for (target in targets) 
        {
            if (Length(target) > Length(register)) 
            {
                fail "Too many targets!";
            }
            let opTargets = Subarray(target, register);
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
        for (target in targets) 
        {
            if (Length(target) > Length(register)) 
            {
                fail "Too many targets!";
            }
            let opTargets = Subarray(target, register);
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
        for (target in targets) 
        {
            if (Length(target) > Length(register)) 
            {
                fail "Too many targets!";
            }
            let opTargets = Subarray(target, register);
            op(opTargets);
        }
    }

}
