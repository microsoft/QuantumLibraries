// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.
namespace Microsoft.Quantum.Synthesis {
    open Microsoft.Quantum.Arithmetic;
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Convert;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Logical;
    open Microsoft.Quantum.Math;

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
        return Where(EqualB(true, _), IntAsBoolArray(value, length));
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
    internal function TBSMain (perm : Int[]) : MCMTMask[] {
        let xs = RangeAsIntArray(0..Length(perm) - 1);
        let gates = new MCMTMask[0];
        return Reversed(Snd(Fold(TBSStep, (perm, gates), xs)));
    }


    /// # Summary
    /// Transform mask of control and target bits to a pair of control qubits and target qubits
    internal function MaskToQubitsPair (qubits : Qubit[], mask : MCMTMask) : (Qubit[], Qubit[]) {
        let n = Length(qubits);
        let controlBits = IntegerBits(mask::ControlMask, n);
        let targetBits = IntegerBits(mask::TargetMask, n);
        let cQubits = Subarray(controlBits, qubits);
        let tQubits = Subarray(targetBits, qubits);

        return (cQubits, tQubits);
    }


    ////////////////////////////////////////////////////////////
    // Public operation                                       //
    ////////////////////////////////////////////////////////////

    /// # Summary
    /// Transformation-based synthesis algorithm.
    ///
    /// This procedure implements the unidirectional transformation based
    /// synthesis approach.  Input is a permutation π over 2ⁿ elements {0, ...,
    /// 2ⁿ-1}, which represents an 𝑛-variable reversible Boolean function.  The
    /// algorithm performs iteratively the following steps:
    ///
    /// 1. Find smallest 𝑥 such that 𝑥 ≠ π(𝑥) = 𝑦
    /// 2. Find multiple-controlled Toffoli gates, which applied to the outputs
    ///    make π(𝑥) = 𝑥 and do not change π(𝑥') for all 𝑥' < 𝑥
    ///
    /// # Input
    /// ## perm
    /// A permutation of 2ⁿ elements starting from 0.
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
        // Translate MCT masks into multiple-controlled multiple-target Toffoli gates.
        let gates = Mapped(MaskToQubitsPair(qubits, _), TBSMain(perm));

        for (gate in gates) {
            let (controls, target) = gate;
            let MultiX = ApplyToEachCA(X, _);
            Controlled MultiX(controls, target);
        }
    } 
}
