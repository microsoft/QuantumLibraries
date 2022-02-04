// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Synthesis {
    open Microsoft.Quantum.Arithmetic;
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Convert;
    open Microsoft.Quantum.Diagnostics;
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
    /// ```qsharp
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
               | [];
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
    internal function UpdatedOutputPattern (pattern : Int, gateMask : MCMTMask) : Int {
        return (pattern &&& gateMask::ControlMask) == gateMask::ControlMask
               ? pattern ^^^ gateMask::TargetMask
               | pattern;
    }


    /// # Summary
    /// Update permutation based according to gate mask.
    internal function UpdatedPermutation (perm : Int[], gateMask : MCMTMask) : Int[] {
        return Mapped(UpdatedOutputPattern(_, gateMask), perm);
    }


    /// # Summary
    /// Computes gate masks to transform perm[x] to x and updates the current
    /// permutation.
    internal function TBSStep (state : (Int[], MCMTMask[]), x : Int) : (Int[], MCMTMask[]) {
        let (perm, gates) = state;
        let y = perm[x];
        let masks = GateMasksForAssignment(x, y);
        let new_perm = Fold(UpdatedPermutation, perm, masks);
        return (new_perm, gates + masks);
    }


    /// # Summary
    /// Compute gate masks to synthesize permutation.
    internal function TBSMain (perm : Int[]) : MCMTMask[] {
        let xs = RangeAsIntArray(0..Length(perm) - 1);
        let gates = [];
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
    /// Permutes the amplitudes in a quantum state given a permutation
    /// using transformation-based synthesis.
    ///
    /// # Description
    /// This procedure implements the unidirectional transformation based
    /// synthesis approach.  Input is a permutation $\pi$ over $2^n$ elements
    /// $\{0, \dots, 2^n-1\}$, which represents an $n$-variable reversible Boolean function. 
    /// The algorithm performs iteratively the following steps:
    ///
    /// 1. Find smallest $x$ such that $x \ne \pi(x) = y$.
    /// 2. Find multiple-controlled Toffoli operations, which applied to the outputs
    ///    make $\pi(x) = x$ and do not change $\pi(x')$ for all $x' < x$
    ///
    /// # Input
    /// ## perm
    /// A permutation of $2^n$ elements starting from 0.
    /// ## qubits
    /// A list of $n$ qubits to which the permutation is applied to.
    ///
    /// # Example
    /// To synthesize a `SWAP` operation:
    /// ```qsharp
    /// using (qubits = Qubit[2]) {
    ///   ApplyPermutationUsingTransformation([0, 2, 1, 3], LittleEndian(qubits));
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
    ///
    /// # See Also
    /// - Microsoft.Quantum.Synthesis.ApplyPermutationUsingDecomposition
    operation ApplyPermutationUsingTransformation(perm : Int[], qubits : LittleEndian) : Unit is Adj + Ctl {
        // Translate MCT masks into multiple-controlled multiple-target Toffoli gates.
        let gates = Mapped(MaskToQubitsPair(qubits!, _), TBSMain(perm));

        Fact(IsPermutation(perm), "perm must be a permutation");
        EqualityFactI(Length(perm), 2^Length(qubits!), $"Length of perm must be {2^Length(qubits!)}");

        for gate in gates {
            let (controls, target) = gate;
            let MultiX = ApplyToEachCA(X, _);
            Controlled MultiX(controls, target);
        }
    }
}
