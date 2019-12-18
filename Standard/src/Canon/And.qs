// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Convert;
    open Microsoft.Quantum.Diagnostics;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Math;
    open Microsoft.Quantum.Measurement;

    /// # Summary
    /// Inverts a given target qubit if and only if both control qubits are in the 1 state,
    /// using measurement to perform the adjoint operation.
    ///
    /// # Description
    /// Inverts `target` if and only if both controls are 1, but assumes that
    /// `target` is in state 0.  The operation has T-count 4, T-depth 2 and
    /// requires no helper qubit, and may therefore be preferable to a CCNOT
    /// operation, if `target` is known to be 0.  The adjoint of this operation
    /// is measurement based and requires no T gates.
    ///
    /// The controlled application of this operation requires no helper qubit,
    /// `2^c` `Rz` gates and is not optimized for depth, where `c` is the number
    /// of overall control qubits including the two controls from the `ApplyAnd`
    /// operation.  The adjoint controlled application requires `2^c - 1` `Rz`
    /// gates (with an angle twice the size of the non-adjoint operation), no
    /// helper qubit and is not optimized for depth.
    ///
    /// # Input
    /// ## control1
    /// First control qubit
    /// ## control2
    /// Second control qubit
    /// ## target
    /// Target auxillary qubit; must be in state 0
    ///
    /// # References
    /// - Cody Jones: "Novel constructions for the fault-tolerant Toffoli gate",
    ///   Phys. Rev. A 87, 022328, 2013
    ///   [arXiv:1212.5069](https://arxiv.org/abs/1212.5069)
    ///   doi:10.1103/PhysRevA.87.022328
    /// - Craig Gidney: "Halving the cost of quantum addition", Quantum 2, page
    ///   74, 2018
    ///   [arXiv:1709.06648](https://arxiv.org/abs/1709.06648)
    ///   doi:10.1103/PhysRevA.85.044302
    /// - Mathias Soeken: "Quantum Oracle Circuits and the Christmas Tree Pattern",
    ///   [Blog article from Decemer 19, 2019](https://msoeken.github.io/blog_qac.html)
    ///   (note: explains the multiple controlled construction)
    operation ApplyAnd(control1 : Qubit, control2 : Qubit, target : Qubit) : Unit {
        body (...) {
            AssertAllZero([target]);
            H(target);
            T(target);
            CNOT(control1, target);
            CNOT(control2, target);
            within {
                CNOT(target, control1);
                CNOT(target, control2);
            }
            apply {
                Adjoint T(control1);
                Adjoint T(control2);
                T(target);
            }
            HY(target);
        }
        adjoint (...) {
            H(target);
            AssertProb([PauliZ], [target], One, 0.5, "Probability of the measurement must be 0.5", 1e-10);
            if (IsResultOne(MResetZ(target))) {
                CZ(control1, control2);
            }
        }
        controlled (controls, ...) {
            _ApplyMultipleControlledAnd(controls + [control1, control2], target);
        }
        adjoint controlled (controls, ...) {
            Adjoint _ApplyMultipleControlledAnd(controls + [control1, control2], target);
        }
    }

    /// # Summary
    /// Inverts a given target qubit if and only if both control qubits are in
    /// the 1 state, with T-depth 1, using measurement to perform the adjoint
    /// operation.
    ///
    /// # Description
    /// Inverts `target` if and only if both controls are 1, but assumes that
    /// `target` is in state 0.  The operation has T-count 4, T-depth 1 and
    /// requires one helper qubit, and may therefore be preferable to a CCNOT
    /// operation, if `target` is known to be 0.  The adjoint of this operation
    /// is measurement based and requires no T gates, and no helper qubit.
    ///
    /// # Input
    /// ## control1
    /// First control qubit
    /// ## control2
    /// Second control qubit
    /// ## target
    /// Target auxillary qubit; must be in state 0
    ///
    /// # References
    /// - Cody Jones: "Novel constructions for the fault-tolerant Toffoli gate",
    ///   Phys. Rev. A 87, 022328, 2013
    ///   [arXiv:1212.5069](https://arxiv.org/abs/1212.5069)
    ///   doi:10.1103/PhysRevA.87.022328
    /// - Peter Selinger: "Quantum circuits of T-depth one",
    ///   Phys. Rev. A 87, 042302, 2013
    ///   [arXiv:1210.0974](https://arxiv.org/abs/1210.0974)
    ///   doi:10.1103/PhysRevA.87.042302
    operation ApplyLowDepthAnd(control1 : Qubit, control2 : Qubit, target : Qubit) : Unit {
        body (...) {
            using (helper = Qubit()) {
                AssertAllZero([target]);
                H(target);
                within {
                    CNOT(target, control1);
                    CNOT(control1, helper);
                    CNOT(control2, helper);
                    CNOT(target, control2);
                }
                apply {
                    Adjoint T(control1);
                    Adjoint T(control2);
                    T(target);
                    T(helper);
                }
                HY(target);
            }
        }
        adjoint (...) {
            Adjoint ApplyAnd(control1, control2, target);
        }
        controlled (controls, ...) {
            _ApplyMultipleControlledLowDepthAnd(controls + [control1, control2], target);
        }
        adjoint controlled (controls, ...) {
            Adjoint _ApplyMultipleControlledLowDepthAnd(controls + [control1, control2], target);
        }
    }

    /// # Summary
    /// Creates Gray code sequences
    ///
    /// # Input
    /// ## n
    /// Number of bits
    ///
    /// # Output
    /// Array of tuples. First value in tuple is value in GrayCode sequence
    /// Second value in tuple is position to change in current value to get
    /// next one.
    ///
    /// # Example
    /// ```Q#
    /// _GrayCode(2); // [(0, 0),(1, 1),(3, 0),(2, 1)]
    /// ```
    function _GrayCode(n : Int) : (Int, Int)[] {
        let N = 1 <<< n;

        mutable res = new (Int, Int)[N];
        mutable j = 0;
        mutable current = IntAsBoolArray(0, n);

        for (i in 0..N - 1) {
            if (i % 2 == 0) {
                set j = 0;
            } else {
                let e = Zip(current, RangeAsIntArray(0..N - 1));
                set j = Snd(Head(Filtered(Fst<Bool, Int>, e))) + 1;
            }

            set j = MaxI(0, Min([j, n - 1]));
            set res w/= i <- (BoolArrayAsInt(current), j);
            if (j < n) {
                set current w/= j <- not current[j];
            }
        }

        return res;
    }

    /// # Summary
    /// Computes the Hamming weight of an integer, i.e., the number of 1s in its
    /// binary expansion.
    ///
    /// # Input
    /// ## number
    /// Number to compute Hamming weight
    /// # Output
    /// Hamming weight of the number
    function _HammingWeightI(number : Int) : Int {
        mutable cnt = number;
        set cnt = (cnt &&& 0x5555555555555555) + ((cnt >>>  1) &&& 0x5555555555555555);
        set cnt = (cnt &&& 0x3333333333333333) + ((cnt >>>  2) &&& 0x3333333333333333);
        set cnt = (cnt &&& 0x0f0f0f0f0f0f0f0f) + ((cnt >>>  4) &&& 0x0f0f0f0f0f0f0f0f);
        set cnt = (cnt &&& 0x00ff00ff00ff00ff) + ((cnt >>>  8) &&& 0x00ff00ff00ff00ff);
        set cnt = (cnt &&& 0x0000ffff0000ffff) + ((cnt >>> 16) &&& 0x0000ffff0000ffff);
        set cnt = (cnt &&& 0x00000000ffffffff) + ((cnt >>> 32) &&& 0x00000000ffffffff);
        return cnt;
    }

    /// # Summary
    /// Returns 1, if `index` has an odd number of 1s and -1, if `index` has an
    /// even number of 1s.
    ///
    /// # Description
    /// Value corresponds to the sign of the coefficient of the Rademacher-Walsh
    /// spectrum of the n-variable AND function for a given assignment that
    /// decides the angle of the rotation.
    ///
    /// # Input
    /// ## index
    /// Input assignment as integer (from 0 to 2^n - 1)
    function _Angle(index : Int) : Int {
        return _HammingWeightI(index) % 2 == 1 ? 1 | -1;
    }

    /// # Summary
    /// Implements a multiple-controlled Toffoli gate, assuming that target
    /// qubit is initialized 0.  The adjoint operation assumes that the target
    /// qubit will be released to 0.
    ///
    /// # Input
    /// ## controls
    /// Control qubits
    /// ## target
    /// Target qubit
    operation _ApplyMultipleControlledAnd(controls : Qubit[], target : Qubit) : Unit {
        body (...) {
            let vars = Length(controls);

            AssertAllZero([target]);

            H(target);

            let code = _GrayCode(vars);
            for (j in 0..Length(code) - 1) {
                let (offset, ctrl) = code[j];
                RFrac(PauliZ, _Angle(offset), vars + 1, target);
                CNOT(controls[ctrl], target);
            }

            HY(target);
        }
        adjoint (...) {
            let vars = Length(controls);

            H(target);
            AssertProb([PauliZ], [target], One, 0.5, "Probability of the measurement must be 0.5", 1e-10);
            if (IsResultOne(MResetZ(target))) {
                for (i in 0..vars - 1) {
                    let start = 1 <<< i;
                    let code = _GrayCode(i);
                    for (j in 0..Length(code) - 1) {
                        let (offset, ctrl) = code[j];
                        RFrac(PauliZ, -_Angle(start + offset), vars, controls[i]);
                        if (i != 0) {
                            CNOT(controls[ctrl], controls[i]);
                        }
                    }
                }
            }
        }
    }

    /// # Summary
    /// Arrange control, target, and helper qubits according to an index
    ///
    /// # Description
    /// Returns a Qubit array with target at index 0, and control i at index
    /// 2^i.  The helper qubits are inserted to all other positions in the
    /// array.
    function _ArrangeQubits(controls : Qubit[], target : Qubit, helper : Qubit[]) : Qubit[] {
        let numControls = Length(controls);
        mutable qs = new Qubit[2^numControls];
        set qs w/= 0 <- target;
        mutable cntC = 0;
        mutable cntH = 0;
        for (i in 1..2^numControls - 1) {
            if (i == (i &&& -i)) {
                set qs w/= i <- controls[cntC];
                set cntC += 1;
            } else {
                set qs w/= i <- helper[cntH];
                set cntH += 1;
            }
        }
        return qs;
    }

    /// # Summary
    /// Implements a multiple-controlled Toffoli gate, assuming that target
    /// qubit is initialized 0.  The adjoint operation assumes that the target
    /// qubit will be released to 0.  Requires a Rz depth of 1, while the number
    /// of helper qubits are exponential in the number of qubits.
    ///
    /// # Input
    /// ## controls
    /// Control qubits
    /// ## target
    /// Target qubit
    operation _ApplyMultipleControlledLowDepthAnd(controls : Qubit[], target : Qubit) : Unit {
        body (...) {
            let vars = Length(controls);
            using (helper = Qubit[2^vars - vars - 1]) {
                let qs = _ArrangeQubits(controls, target, helper);

                AssertAllZero([target]);
                H(target);

                within {
                    // initialize helper lines with control lines based on LSB
                    for (i in 3..2^vars - 1) {
                        let lsb = i &&& -i;
                        if (i != lsb) { // i is power of 2
                            CNOT(qs[lsb], qs[i]);
                        }
                    }
                    // target to control
                    ApplyToEachA(CNOT(target, _), controls);
                    // copy remainder (without LSB)
                    for (i in 3..2^vars - 1) {
                        let lsb = i &&& -i;
                        if (i != lsb) {
                            CNOT(qs[i - lsb], qs[i]);
                        }
                    }
                } apply {
                    for (i in IndexRange(qs)) {
                        RFrac(PauliZ, _Angle(i), vars + 1, qs[i]);
                    }
                }

                HY(target);
            }
        }
        adjoint (...) {
            let vars = Length(controls);

            H(target);
            AssertProb([PauliZ], [target], One, 0.5, "Probability of the measurement must be 0.5", 1e-10);
            if (IsResultOne(MResetZ(target))) {
                using (helper = Qubit[2^vars - vars - 1]) {
                    let qs = _ArrangeQubits(controls, target, helper);
                    within {
                        // this is a bit easier than in the compute part, since
                        // the target qubit does not have to be copied over to
                        // the control lines.  Therefore, the two LSB CNOT parts
                        // can be merged into a single loop.
                        for (i in 3..2^vars - 1) {
                            let lsb = i &&& -i;
                            if (i != lsb) {
                                CNOT(qs[lsb], qs[i]);
                                CNOT(qs[i - lsb], qs[i]);
                            }
                        }
                    } apply {
                        for (i in 1..2^vars - 1) {
                            RFrac(PauliZ, -_Angle(i), vars, qs[i]);
                        }
                    }
                }
            }
        }
    }
}
