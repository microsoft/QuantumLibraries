// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {

    /// # Summary
    /// Applies an operation conditioned on a classical bit.
    ///
    /// # Description
    /// Given an operation `op` and a bit value `bit`, applies `op` to the `target`
    /// if `bit` is `true`. If `false`, nothing happens to the `target`.
    ///
    /// # Input
    /// ## op
    /// An operation to be conditionally applied.
    /// ## bit
    /// a boolean that controls whether op is applied or not.
    /// ## target
    /// The input to which the operation is applied.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The input type of the operation to be conditionally applied.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApplyIfC
    /// - Microsoft.Quantum.Canon.ApplyIfA
    /// - Microsoft.Quantum.Canon.ApplyIfCA
    ///
    /// # Example
    /// The following prepares a register of qubits into a computational basis
    /// state represented by a classical bit string given as an array of `Bool`
    /// values:
    /// ```Q#
    /// let bitstring = [true, false, true];
    /// using (register = Qubit(3)) {
    ///     ApplyToEach(ApplyIf(X, _, _), Zip(bitstring, register));
    ///     // register should now be in the state |101⟩.
    ///     ...
    /// }
    /// ```
    operation ApplyIf<'T> (op : ('T => Unit), bit : Bool, target : 'T) : Unit {
        if (bit) {
            op(target);
        }
    }

    /// # Summary
    /// Applies a controllable operation conditioned on a classical bit.
    ///
    /// # Description
    /// Given an operation `op` and a bit value `bit`, applies `op` to the `target`
    /// if `bit` is `true`. If `false`, nothing happens to the `target`.
    /// The suffix `C` indicates that the operation to be applied is controllable.
    ///
    /// # Input
    /// ## op
    /// An operation to be conditionally applied.
    /// ## bit
    /// a boolean that controls whether op is applied or not.
    /// ## target
    /// The input to which the operation is applied.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The input type of the operation to be conditionally applied.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApplyIfC
    /// - Microsoft.Quantum.Canon.ApplyIfA
    /// - Microsoft.Quantum.Canon.ApplyIfCA
    ///
    /// # Example
    /// The following prepares a register of qubits into a computational basis
    /// state represented by a classical bit string given as an array of `Bool`
    /// values:
    /// ```Q#
    /// let bitstring = [true, false, true];
    /// using (register = Qubit(3)) {
    ///     ApplyToEach(ApplyIf(X, _, _), Zip(bitstring, register));
    ///     // register should now be in the state |101⟩.
    ///     ...
    /// }
    /// ```
    operation ApplyIfC<'T> (op : ('T => Unit is Ctl), bit : Bool, target : 'T) : Unit is Ctl {
        if (bit) {
            op(target);
        }
    }

    /// # Summary
    /// Applies a adjointable operation conditioned on a classical bit.
    ///
    /// # Description
    /// Given an operation `op` and a bit value `bit`, applies `op` to the `target`
    /// if `bit` is `true`. If `false`, nothing happens to the `target`.
    /// The suffix `A` indicates that the operation to be applied is adjointable.
    ///
    /// # Input
    /// ## op
    /// An operation to be conditionally applied.
    /// ## bit
    /// a boolean that controls whether op is applied or not.
    /// ## target
    /// The input to which the operation is applied.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The input type of the operation to be conditionally applied.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApplyIfC
    /// - Microsoft.Quantum.Canon.ApplyIfA
    /// - Microsoft.Quantum.Canon.ApplyIfCA
    ///
    /// # Example
    /// The following prepares a register of qubits into a computational basis
    /// state represented by a classical bit string given as an array of `Bool`
    /// values:
    /// ```Q#
    /// let bitstring = [true, false, true];
    /// using (register = Qubit(3)) {
    ///     ApplyToEach(ApplyIf(X, _, _), Zip(bitstring, register));
    ///     // register should now be in the state |101⟩.
    ///     ...
    /// }
    /// ```
    operation ApplyIfA<'T> (op : ('T => Unit is Adj), bit : Bool, target : 'T) : Unit is Adj {
        if (bit) {
            op(target);
        }
    }

    /// # Summary
    /// Applies a unitary operation conditioned on a classical bit.
    ///
    /// # Description
    /// Given an operation `op` and a bit value `bit`, applies `op` to the `target`
    /// if `bit` is `true`. If `false`, nothing happens to the `target`.
    /// The suffix `CA` indicates that the operation to be applied is unitary
    /// (controllable and adjointable).
    ///
    /// # Input
    /// ## op
    /// An operation to be conditionally applied.
    /// ## bit
    /// a boolean that controls whether op is applied or not.
    /// ## target
    /// The input to which the operation is applied.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The input type of the operation to be conditionally applied.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApplyIfC
    /// - Microsoft.Quantum.Canon.ApplyIfA
    /// - Microsoft.Quantum.Canon.ApplyIfCA
    ///
    /// # Example
    /// The following prepares a register of qubits into a computational basis
    /// state represented by a classical bit string given as an array of `Bool`
    /// values:
    /// ```Q#
    /// let bitstring = [true, false, true];
    /// using (register = Qubit(3)) {
    ///     ApplyToEach(ApplyIf(X, _, _), Zip(bitstring, register));
    ///     // register should now be in the state |101⟩.
    ///     ...
    /// }
    /// ```
    operation ApplyIfCA<'T> (op : ('T => Unit is Ctl + Adj), bit : Bool, target : 'T) : Unit is Ctl + Adj {
        if (bit) {
            op(target);
        }
    }

    /// # Summary
    /// Applies an operation conditioned on a classical result value being one.
    ///
    /// # Description
    /// Given an operation `op` and a result value `result`, applies `op` to the `target`
    /// if `result` is `One`. If `Zero`, nothing happens to the `target`.
    ///
    /// # Input
    /// ## op
    /// An operation to be conditionally applied.
    /// ## result
    /// A measurement result that controls whether op is applied or not.
    /// ## target
    /// The input to which the operation is applied.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The input type of the operation to be conditionally applied.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApplyIfOneC
    /// - Microsoft.Quantum.Canon.ApplyIfOneA
    /// - Microsoft.Quantum.Canon.ApplyIfOneCA
    operation ApplyIfOne<'T> (result : Result, (op : ('T => Unit), target : 'T)) : Unit {
        if (result == One) {
            op(target);
        }
    }

    /// # Summary
    /// Applies a controllable operation conditioned on a classical result value being one.
    ///
    /// # Description
    /// Given an operation `op` and a result value `result`, applies `op` to the `target`
    /// if `result` is `One`. If `Zero`, nothing happens to the `target`.
    /// The suffix `C` indicates that the operation to be applied is controllable.
    ///
    /// # Input
    /// ## op
    /// An operation to be conditionally applied.
    /// ## result
    /// A measurement result that controls whether op is applied or not.
    /// ## target
    /// The input to which the operation is applied.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The input type of the operation to be conditionally applied.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApplyIfOneC
    /// - Microsoft.Quantum.Canon.ApplyIfOneA
    /// - Microsoft.Quantum.Canon.ApplyIfOneCA
    operation ApplyIfOneC<'T> (result : Result, (op : ('T => Unit is Ctl), target : 'T)) : Unit is Ctl {
        if (result == One) {
            op(target);
        }
    }

    /// # Summary
    /// Applies an adjointable operation conditioned on a classical result value being one.
    ///
    /// # Description
    /// Given an operation `op` and a result value `result`, applies `op` to the `target`
    /// if `result` is `One`. If `Zero`, nothing happens to the `target`.
    /// The suffix `A` indicates that the operation to be applied is adjointable.
    ///
    /// # Input
    /// ## op
    /// An operation to be conditionally applied.
    /// ## result
    /// A measurement result that controls whether op is applied or not.
    /// ## target
    /// The input to which the operation is applied.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The input type of the operation to be conditionally applied.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApplyIfOneC
    /// - Microsoft.Quantum.Canon.ApplyIfOneA
    /// - Microsoft.Quantum.Canon.ApplyIfOneCA
    operation ApplyIfOneA<'T> (result : Result, (op : ('T => Unit is Adj), target : 'T)) : Unit is Adj {
        if (result == One) {
            op(target);
        }
    }

    /// # Summary
    /// Applies a unitary operation conditioned on a classical result value being one.
    ///
    /// # Description
    /// Given an operation `op` and a result value `result`, applies `op` to the `target`
    /// if `result` is `One`. If `Zero`, nothing happens to the `target`.
    /// The suffix `CA` indicates that the operation to be applied is unitary
    /// (controllable and adjointable).
    ///
    /// # Input
    /// ## op
    /// An operation to be conditionally applied.
    /// ## result
    /// A measurement result that controls whether op is applied or not.
    /// ## target
    /// The input to which the operation is applied.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The input type of the operation to be conditionally applied.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApplyIfOneC
    /// - Microsoft.Quantum.Canon.ApplyIfOneA
    /// - Microsoft.Quantum.Canon.ApplyIfOneCA
    operation ApplyIfOneCA<'T> (result : Result, (op : ('T => Unit is Adj + Ctl), target : 'T)) : Unit is Adj + Ctl {
        if (result == One) {
            op(target);
        }
    }

    /// # Summary
    /// Applies an operation conditioned on a classical result value being zero.
    ///
    /// # Description
    /// Given an operation `op` and a result value `result`, applies `op` to the `target`
    /// if `result` is `Zero`. If `One`, nothing happens to the `target`.
    ///
    /// # Input
    /// ## op
    /// An operation to be conditionally applied.
    /// ## result
    /// A measurement result that controls whether op is applied or not.
    /// ## target
    /// The input to which the operation is applied.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The input type of the operation to be conditionally applied.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApplyIfZeroC
    /// - Microsoft.Quantum.Canon.ApplyIfZeroA
    /// - Microsoft.Quantum.Canon.ApplyIfZeroCA
    operation ApplyIfZero<'T> (result : Result, (op : ('T => Unit), target : 'T)) : Unit {
        if (result == Zero) {
            op(target);
        }
    }

    /// # Summary
    /// Applies a controllable operation conditioned on a classical result value being zero.
    ///
    /// # Description
    /// Given an operation `op` and a result value `result`, applies `op` to the `target`
    /// if `result` is `Zero`. If `One`, nothing happens to the `target`.
    /// The suffix `C` indicates that the operation to be applied is controllable.
    ///
    /// # Input
    /// ## op
    /// An operation to be conditionally applied.
    /// ## result
    /// A measurement result that controls whether op is applied or not.
    /// ## target
    /// The input to which the operation is applied.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The input type of the operation to be conditionally applied.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApplyIfZeroC
    /// - Microsoft.Quantum.Canon.ApplyIfZeroA
    /// - Microsoft.Quantum.Canon.ApplyIfZeroCA
    operation ApplyIfZeroC<'T> (result : Result, (op : ('T => Unit is Ctl), target : 'T)) : Unit is Ctl {
        if (result == Zero) {
            op(target);
        }
    }

    /// # Summary
    /// Applies an adjointable operation conditioned on a classical result value being zero.
    ///
    /// # Description
    /// Given an operation `op` and a result value `result`, applies `op` to the `target`
    /// if `result` is `Zero`. If `One`, nothing happens to the `target`.
    /// The suffix `A` indicates that the operation to be applied is adjointable.
    ///
    /// # Input
    /// ## op
    /// An operation to be conditionally applied.
    /// ## result
    /// A measurement result that controls whether op is applied or not.
    /// ## target
    /// The input to which the operation is applied.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The input type of the operation to be conditionally applied.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApplyIfZeroC
    /// - Microsoft.Quantum.Canon.ApplyIfZeroA
    /// - Microsoft.Quantum.Canon.ApplyIfZeroCA
    operation ApplyIfZeroA<'T> (result : Result, (op : ('T => Unit is Adj), target : 'T)) : Unit is Adj {
        if (result == Zero) {
            op(target);
        }
    }

    /// # Summary
    /// Applies a unitary operation conditioned on a classical result value being zero.
    ///
    /// # Description
    /// Given an operation `op` and a result value `result`, applies `op` to the `target`
    /// if `result` is `Zero`. If `One`, nothing happens to the `target`.
    /// The suffix `CA` indicates that the operation to be applied is unitary
    /// (controllable and adjointable).
    ///
    /// # Input
    /// ## op
    /// An operation to be conditionally applied.
    /// ## result
    /// A measurement result that controls whether op is applied or not.
    /// ## target
    /// The input to which the operation is applied.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The input type of the operation to be conditionally applied.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApplyIfZeroC
    /// - Microsoft.Quantum.Canon.ApplyIfZeroA
    /// - Microsoft.Quantum.Canon.ApplyIfZeroCA
    operation ApplyIfZeroCA<'T> (result : Result, (op : ('T => Unit is Adj + Ctl), target : 'T)) : Unit is Adj + Ctl {
        if (result == Zero) {
            op(target);
        }
    }

    /// # Summary
    /// Applies one of two operations, depending on the value of a classical
    /// result.
    ///
    /// # Description
    /// Given a result `result`, applies the operation `zeroOp` with `zeroInput` as
    /// its input when `result` is equal to `Zero`, and applies `oneOp(oneInput)`
    /// when `result == One`.
    ///
    /// # Input
    /// ## result
    /// The measurement result used to determine if `zeroOp` or `oneOp` is
    /// applied.
    /// ## zeroOp
    /// The operation to be applied when `result == Zero`.
    /// ## zeroInput
    /// The input to be provided to `zeroOp` when `result == Zero`.
    /// ## oneOp
    /// The operation to be applied when `result == One`.
    /// ## oneInput
    /// The input to be provided to `oneOp` when `result == One`.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The input type of the operation `zeroOp` to be conditionally applied.
    /// ## 'U
    /// The input type of the operation `oneOp` to be conditionally applied.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApplyIfZero
    /// - Microsoft.Quantum.Canon.ApplyIfOne
    /// - Microsoft.Quantum.Canon.ApplyIfElseRC
    /// - Microsoft.Quantum.Canon.ApplyIfElseRA
    /// - Microsoft.Quantum.Canon.ApplyIfElseRCA
    operation ApplyIfElseR<'T,'U> (result : Result, (zeroOp : ('T => Unit), zeroInput : 'T), (oneOp : ('U => Unit), oneInput : 'U)) : Unit {
        if (result == Zero) {
            zeroOp(zeroInput);
        } else {
            oneOp(oneInput);
        }
    }

    /// # Summary
    /// Applies one of two controllable operations, depending on the value of a
    /// classical result.
    ///
    /// # Description
    /// Given a result `result`, applies the operation `zeroOp` with `zeroInput` as
    /// its input when `result` is equal to `Zero`, and applies `oneOp(oneInput)`
    /// when `result == One`.
    ///
    /// # Input
    /// ## result
    /// The measurement result used to determine if `zeroOp` or `oneOp` is
    /// applied.
    /// ## zeroOp
    /// The controllable operation to be applied when `result == Zero`.
    /// ## zeroInput
    /// The input to be provided to `zeroOp` when `result == Zero`.
    /// ## oneOp
    /// The controllable operation to be applied when `result == One`.
    /// ## oneInput
    /// The input to be provided to `oneOp` when `result == One`.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The input type of the operation `zeroOp` to be conditionally applied.
    /// ## 'U
    /// The input type of the operation `oneOp` to be conditionally applied.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApplyIfZero
    /// - Microsoft.Quantum.Canon.ApplyIfOne
    /// - Microsoft.Quantum.Canon.ApplyIfElseRC
    /// - Microsoft.Quantum.Canon.ApplyIfElseRA
    /// - Microsoft.Quantum.Canon.ApplyIfElseRCA
    operation ApplyIfElseRC<'T,'U> (result : Result, (zeroOp : ('T => Unit is Ctl), zeroInput : 'T), (oneOp : ('U => Unit is Ctl), oneInput : 'U)) : Unit is Ctl {
        if (result == Zero) {
            zeroOp(zeroInput);
        } else {
            oneOp(oneInput);
        }
    }

    /// # Summary
    /// Applies one of two adjointable operations, depending on the value of a
    /// classical result.
    ///
    /// # Description
    /// Given a result `result`, applies the operation `zeroOp` with `zeroInput` as
    /// its input when `result` is equal to `Zero`, and applies `oneOp(oneInput)`
    /// when `result == One`.
    ///
    /// # Input
    /// ## result
    /// The measurement result used to determine if `zeroOp` or `oneOp` is
    /// applied.
    /// ## zeroOp
    /// The adjointable operation to be applied when `result == Zero`.
    /// ## zeroInput
    /// The input to be provided to `zeroOp` when `result == Zero`.
    /// ## oneOp
    /// The adjointable operation to be applied when `result == One`.
    /// ## oneInput
    /// The input to be provided to `oneOp` when `result == One`.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The input type of the operation `zeroOp` to be conditionally applied.
    /// ## 'U
    /// The input type of the operation `oneOp` to be conditionally applied.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApplyIfZero
    /// - Microsoft.Quantum.Canon.ApplyIfOne
    /// - Microsoft.Quantum.Canon.ApplyIfElseRC
    /// - Microsoft.Quantum.Canon.ApplyIfElseRA
    /// - Microsoft.Quantum.Canon.ApplyIfElseRCA
    operation ApplyIfElseRA<'T,'U> (result : Result, (zeroOp : ('T => Unit is Adj), zeroInput : 'T), (oneOp : ('U => Unit is Adj), oneInput : 'U)) : Unit is Adj {
        if (result == Zero) {
            zeroOp(zeroInput);
        } else {
            oneOp(oneInput);
        }
    }

    /// # Summary
    /// Applies one of two unitary operations, depending on the value of a
    /// classical result.
    ///
    /// # Description
    /// Given a result `result`, applies the operation `zeroOp` with `zeroInput` as
    /// its input when `result` is equal to `Zero`, and applies `oneOp(oneInput)`
    /// when `result == One`.
    ///
    /// # Input
    /// ## result
    /// The measurement result used to determine if `zeroOp` or `oneOp` is
    /// applied.
    /// ## zeroOp
    /// The unitary operation to be applied when `result == Zero`.
    /// ## zeroInput
    /// The input to be provided to `zeroOp` when `result == Zero`.
    /// ## oneOp
    /// The unitary operation to be applied when `result == One`.
    /// ## oneInput
    /// The input to be provided to `oneOp` when `result == One`.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The input type of the operation `zeroOp` to be conditionally applied.
    /// ## 'U
    /// The input type of the operation `oneOp` to be conditionally applied.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApplyIfZero
    /// - Microsoft.Quantum.Canon.ApplyIfOne
    /// - Microsoft.Quantum.Canon.ApplyIfElseRC
    /// - Microsoft.Quantum.Canon.ApplyIfElseRA
    /// - Microsoft.Quantum.Canon.ApplyIfElseRCA
    operation ApplyIfElseRCA<'T,'U> (result : Result, (zeroOp : ('T => Unit is Adj + Ctl), zeroInput : 'T), (oneOp : ('U => Unit is Adj + Ctl), oneInput : 'U)) : Unit is Adj + Ctl {
        if (result == Zero) {
            zeroOp(zeroInput);
        } else {
            oneOp(oneInput);
        }
    }

    /// # Summary
    /// Applies one of two operations, depending on the value of a classical
    /// bit.
    ///
    /// # Description
    /// Given a bit `bit`, applies the operation `trueOp` with `trueInput` as
    /// its input when `bit` is `true`, and applies `falseOp(falseInput)`
    /// when `bit` is `false`.
    ///
    /// # Input
    /// ## bit
    /// The boolean value used to determine whether `trueOp` or `falseOp` is
    /// applied.
    /// ## trueOp
    /// The operation to be applied when `bit` is `true`.
    /// ## trueInput
    /// The input to be provided to `trueOp` when `bit` is `true`.
    /// ## falseOp
    /// The operation to be applied when `bit` is `false`.
    /// ## falseInput
    /// The input to be provided to `falseOp` when `bit` is `false`.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The input type of the operation `trueOp` to be conditionally applied.
    /// ## 'U
    /// The input type of the operation `falseOp` to be conditionally applied.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApplyIfZero
    /// - Microsoft.Quantum.Canon.ApplyIfOne
    /// - Microsoft.Quantum.Canon.ApplyIfElseRC
    /// - Microsoft.Quantum.Canon.ApplyIfElseRA
    /// - Microsoft.Quantum.Canon.ApplyIfElseRCA
    operation ApplyIfElseB<'T,'U> (bit : Bool, (trueOp : ('T => Unit), trueInput : 'T), (falseOp : ('U => Unit), falseInput : 'U)) : Unit {
        if (bit) {
            trueOp(trueInput);
        } else {
            falseOp(falseInput);
        }
    }

    /// # Summary
    /// Applies one of two controllable operations, depending on the value of a
    /// classical bit.
    ///
    /// # Description
    /// Given a bit `bit`, applies the operation `trueOp` with `trueInput` as
    /// its input when `bit` is `true`, and applies `falseOp(falseInput)`
    /// when `bit` is `false`.
    ///
    /// # Input
    /// ## bit
    /// The boolean value used to determine whether `trueOp` or `falseOp` is
    /// applied.
    /// ## trueOp
    /// The controllable operation to be applied when `bit` is `true`.
    /// ## trueInput
    /// The input to be provided to `trueOp` when `bit` is `true`.
    /// ## falseOp
    /// The controllable operation to be applied when `bit` is `false`.
    /// ## falseInput
    /// The input to be provided to `falseOp` when `bit` is `false`.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The input type of the operation `trueOp` to be conditionally applied.
    /// ## 'U
    /// The input type of the operation `falseOp` to be conditionally applied.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApplyIfZero
    /// - Microsoft.Quantum.Canon.ApplyIfOne
    /// - Microsoft.Quantum.Canon.ApplyIfElseRC
    /// - Microsoft.Quantum.Canon.ApplyIfElseRA
    /// - Microsoft.Quantum.Canon.ApplyIfElseRCA
    operation ApplyIfElseBC<'T,'U> (bit : Bool, (trueOp : ('T => Unit is Ctl), trueInput : 'T), (falseOp : ('U => Unit is Ctl), falseInput : 'U)) : Unit is Ctl {
        if (bit) {
            trueOp(trueInput);
        } else {
            falseOp(falseInput);
        }
    }

    /// # Summary
    /// Applies one of two adjointable operations, depending on the value of a
    /// classical bit.
    ///
    /// # Description
    /// Given a bit `bit`, applies the operation `trueOp` with `trueInput` as
    /// its input when `bit` is `true`, and applies `falseOp(falseInput)`
    /// when `bit` is `false`.
    ///
    /// # Input
    /// ## bit
    /// The boolean value used to determine if `trueOp` or `falseOp` is
    /// applied.
    /// ## trueOp
    /// The adjointable operation to be applied when `bit` is `true`.
    /// ## trueInput
    /// The input to be provided to `trueOp` when `bit` is `true`.
    /// ## falseOp
    /// The adjointable operation to be applied when `bit` is `false`.
    /// ## falseInput
    /// The input to be provided to `falseOp` when `bit` is `false`.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The input type of the operation `trueOp` to be conditionally applied.
    /// ## 'U
    /// The input type of the operation `falseOp` to be conditionally applied.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApplyIfZero
    /// - Microsoft.Quantum.Canon.ApplyIfOne
    /// - Microsoft.Quantum.Canon.ApplyIfElseRC
    /// - Microsoft.Quantum.Canon.ApplyIfElseRA
    /// - Microsoft.Quantum.Canon.ApplyIfElseRCA
    operation ApplyIfElseBA<'T,'U> (bit : Bool, (trueOp : ('T => Unit is Adj), trueInput : 'T), (falseOp : ('U => Unit is Adj), falseInput : 'U)) : Unit is Adj {
        if (bit) {
            trueOp(trueInput);
        } else {
            falseOp(falseInput);
        }
    }

    /// # Summary
    /// Applies one of two unitary operations, depending on the value of a
    /// classical bit.
    ///
    /// # Description
    /// Given a bit `bit`, applies the operation `trueOp` with `trueInput` as
    /// its input when `bit` is `true`, and applies `falseOp(falseInput)`
    /// when `bit` is `false`.
    ///
    /// # Input
    /// ## bit
    /// The boolean value used to determine if `trueOp` or `falseOp` is
    /// applied.
    /// ## trueOp
    /// The unitary operation to be applied when `bit` is `true`.
    /// ## trueInput
    /// The input to be provided to `trueOp` when `bit` is `true`.
    /// ## falseOp
    /// The unitary operation to be applied when `bit` is `false`.
    /// ## falseInput
    /// The input to be provided to `falseOp` when `bit` is `false`.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The input type of the operation `trueOp` to be conditionally applied.
    /// ## 'U
    /// The input type of the operation `falseOp` to be conditionally applied.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApplyIfZero
    /// - Microsoft.Quantum.Canon.ApplyIfOne
    /// - Microsoft.Quantum.Canon.ApplyIfElseRC
    /// - Microsoft.Quantum.Canon.ApplyIfElseRA
    /// - Microsoft.Quantum.Canon.ApplyIfElseRCA
    operation ApplyIfElseBCA<'T,'U> (bit : Bool, (trueOp : ('T => Unit is Adj + Ctl), trueInput : 'T), (falseOp : ('U => Unit is Adj + Ctl), falseInput : 'U)) : Unit is Adj + Ctl {
        if (bit) {
            trueOp(trueInput);
        } else {
            falseOp(falseInput);
        }
    }

}
