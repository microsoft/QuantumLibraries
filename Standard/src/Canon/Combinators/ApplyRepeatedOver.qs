// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {
    open Microsoft.Quantum.Diagnostics;
    open Microsoft.Quantum.Arrays;
    
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // Helpers to repeatedly apply functions over qubit arrays
    ///////////////////////////////////////////////////////////////////////////////////////////////

    operation ApplySeriesOfOps(listOfOps : (Qubit[] => Unit)[], targets : Int[][], register : Qubit[]) : Unit {
        if (Length(listOfOps) != Length(targets)) {
            fail "The number of ops and number of targets do not match!";
        }
        for (index in 0..Length(listOfOps) - 1) {
            if (Length(targets[index]) > Length(register)) {
                fail "There are too many targets!";
            }
            let opToApply = listOfOps[index];
            let qubitTargets = Subarray(targets[index, register]);
            opToApply(qubitTargets);
        }
    }

    /// # Summary
    /// Applies a multiply controlled version of a singly controlled
    /// operation.
    /// The modifier `C` indicates that the single-qubit operation is controllable.
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

    operation PermuteQubits(targets : Int[][], register : Qubit[]) : Unit {
        ApplyOpRepeatedlyOverCA(ApplyToFirstTwoQubits(SWAP, _), targets, register);
    }

}