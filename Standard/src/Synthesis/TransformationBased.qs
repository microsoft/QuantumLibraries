// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.
namespace Microsoft.Quantum.Synthesis {
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Math;
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Arithmetic;
    open Microsoft.Quantum.Convert;

    ////////////////////////////////////////////////////////////
    // Transformation-based synthesis                         //
    ////////////////////////////////////////////////////////////

    /// # Summary
    /// A type to represent a multiple-controlled multiple-target Toffoli gate.
    ///
    /// The first integer is a bit mask for control lines.  Bit indexes which
    /// are set correspond to control line indexes.
    ///
    /// The second integer is a bit mask for target lines.  Bit indexes which
    /// are set correspond to target line indexes.
    ///
    /// The bit indexes of both integers must be disjoint.
    internal newtype MCMTMask = (
        ControlMask : Int, TargetMask : Int
    );

    /// # Summary
    /// A type to represent a multiple-controlled Toffoli gate.
    ///
    /// The first value is an array of qubits for the control lines, the second
    /// value is a qubit for the target line.
    ///
    /// The target cannot be contained in the control lines.
    ///
    /// # Example
    /// ```Q#
    /// using (qubits = Qubit[4]) {
    ///   let not_gate = MCTGate(new Qubit[0], qubits[0]);
    ///   let cnot_gate = MCTGate([qubits[0]], qubits[1]);
    ///   let toffoli_gate = MCTGate([qubits[0], qubits[1]], qubits[2]);
    ///   let gate = MCTGate([qubits[0], qubits[1], qubits[2]], qubits[3]);
    /// }
    /// ```
    internal newtype MCTGate = (
        Controls : Qubit[],
        Target : Qubit
    );


    // Some helper functions

    /// # Summary
    /// Checks whether bit at position is set in nonnegative number.
    ///
    /// # Input
    /// ## value
    /// A nonnegative number.
    /// ## position
    /// A bit position starting from 0.
    ///
    /// # Output
    /// Returns true, if in the binary expansion of `value` the bit at position
    /// `position` is 1, otherwise false.
    ///
    /// # Example
    /// ```Q#
    /// IsBitSet(23, 0); // true, since 23 is 10111 in binary
    /// IsBitSet(23, 3); // false
    /// ```
    ///
    /// # Note
    /// Implementation with right-shift did not work
    internal function IsBitSet (value : Int, position : Int) : Bool {
        return (value &&& 2 ^ position) == 2 ^ position;
    }


    /// # Summary
    /// Returns all positions in which bits of an integer are set.
    ///
    /// # Input
    /// ## value
    /// A nonnegative number.
    /// ## length
    /// The number of bits in the binary expansion of `value`.
    ///
    /// # Output
    /// An array containing all bit positions (starting from 0) that are 1 in
    /// the binary expansion of `value` considering all bits up to position
    /// `length - 1`.  All positions are ordered in the array by position in an
    /// increasing order.
    ///
    /// # Example
    /// ```Q#
    /// IntegerBits(23, 5); // [0, 1, 2, 4]
    /// IntegerBits(10, 4); // [1, 3]
    /// ```
    internal function IntegerBits (value : Int, length : Int) : Int[] {
        return Filtered(IsBitSet(value, _), RangeAsIntArray(0..length - 1));
    }


    /// # Summary
    /// Constructs a MCMTMask type as a singleton array if targets is not 0,
    /// otherwise returns an empty array.
    internal function GateMask (controls : Int, targets : Int) : MCMTMask[] {
        return targets != 0
               ? [MCMTMask(controls, targets)]
               | new MCMTMask[0];
    }


    /// # Summary
    /// Computes up to two MCMT masks to transform y to x.
    internal function GateMasksForAssignment (x : Int, y : Int) : MCMTMask[] {
        let m01 = x &&& ~~~y;
        let m10 = y &&& ~~~x;
        return GateMask(y, m01) + GateMask(x, m10);
    }


    /// # Summary
    /// Update an output pattern according to gate mask.
    internal function UpdateOutputPattern (pattern : Int, gateMask : MCMTMask) : Int {
        return (pattern &&& gateMask::ControlMask) == gateMask::ControlMask
               ? pattern ^^^ gateMask::TargetMask
               | pattern;
    }


    /// # Summary
    /// Update permutation based according to gate mask.
    internal function UpdatePermutation (perm : Int[], gateMask : MCMTMask) : Int[] {
        return Mapped(UpdateOutputPattern(_, gateMask), perm);
    }


    /// # Summary
    /// Computes gate masks to transform perm[x] to x and updates the current
    /// permutation.
    internal function TBSStep (state : (Int[], MCMTMask[]), x : Int) : (Int[], MCMTMask[]) {
        let (perm, gates) = state;
        let y = perm[x];
        let masks = GateMasksForAssignment(x, y);
        let new_perm = Fold(UpdatePermutation, perm, masks);
        return (new_perm, gates + masks);
    }


    /// # Summary
    /// Compute gate masks to synthesize permutation.
    internal function TBSMain(perm : Int[]) : MCMTMask[] {
        let xs = RangeAsIntArray(0..Length(perm) - 1);
        let gates = new MCMTMask[0];
        return Reversed(Snd(Fold(TBSStep, (perm, gates), xs)));
    }


    /// # Summary
    /// Translate MCT masks into multiple-controlled Toffoli gates (with single
    /// targets).
    internal function GateMasksToToffoliGates (qubits : Qubit[], masks : MCMTMask[]) : MCTGate[] {
        mutable result = new MCTGate[0];
        let n = Length(qubits);

        for (mask in masks) {
            let (controls, targets) = mask!;
            let controlBits = IntegerBits(controls, n);
            let targetBits = IntegerBits(targets, n);
            let cQubits = Subarray(controlBits, qubits);
            let tQubits = Subarray(targetBits, qubits);

            for (t in tQubits) {
                set result += [MCTGate(cQubits, t)];
            }
        }

        return result;
    }

    ////////////////////////////////////////////////////////////
    // Public operation                                       //
    ////////////////////////////////////////////////////////////

    /// # Summary
    /// Transformation-based synthesis algorithm.
    ///
    /// This procedure implements the unidirectional transformation based
    /// synthesis approach.  Input is a permutation π over 2^𝑛 elements {0, ...,
    /// 2^𝑛-1}, which represents an 𝑛-variable reversible Boolean function.  The
    /// algorithm performs iteratively the following steps:
    ///
    /// 1. Find smallest 𝑥 such that 𝑥 ≠ π(𝑥) = 𝑦
    /// 2. Find multiple-controlled Toffoli gates, which applied to the outputs
    ///    make π(𝑥) = 𝑥 and do not change π(𝑥') for all 𝑥' < 𝑥
    ///
    /// # Input
    /// ## perm
    /// A permutation of 2^𝑛 elements starting from 0.
    /// ## qubits
    /// A list of 𝑛 qubits where the Toffoli gates are being applied to.  Note
    /// that the algorithm does not apply the gates.  But only prepares the
    /// Toffoli gates.
    ///
    /// # Example
    /// ```Q#
    /// using (qubits = Qubit[3]) {
    ///   ApplyPermutationTransformationBased([0, 2, 1, 3], qubits); // synthesize SWAP operation
    /// }
    /// ```
    ///
    /// # References
    /// - [*D. Michael Miller*, *Dmitri Maslov*, *Gerhard W. Dueck*,
    ///    Proc. DAC 2003, IEEE, pp. 318-323,
    ///    2003](https://doi.org/10.1145/775832.775915)
    /// - [*Mathias Soeken*, *Gerhard W. Dueck*, *D. Michael Miller*,
    ///    Proc. RC 2016, Springer, pp. 307-321,
    ///    2016](https://doi.org/10.1007/978-3-319-40578-0_22)
    operation ApplyPermutationTransformationBased(perm : Int[], qubits : Qubit[]) : Unit is Adj + Ctl {
        let gates = GateMasksToToffoliGates(qubits, TBSMain(perm));

        for (gate in gates) {
            Controlled X(gate::Controls, gate::Target);
        }
    } 
}
