// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {
    open Microsoft.Quantum.Intrinsic;

    /// # Summary
    /// Given two operations, applies one as conjugated with the other.
    ///
    /// # Description
    /// Given two operations, respectively described by unitary operators $U$
    /// and $V$, applies them in the sequence $U^{\dagger} V U$. That is,
    /// this operation implements the unitary operator given by $V$ conjugated
    /// with $U$.
    ///
    /// # Input
    /// ## outerOperation
    /// The operation $U$ that should be used to conjugate $V$. Note that the
    /// outer operation $U$ needs to be adjointable, but does not
    /// need to be controllable.
    /// ## innerOperation
    /// The operation $V$ being conjugated.
    /// ## target
    /// The input to be provided to the outer and inner operations.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The target on which each of the inner and outer operations act.
    ///
    /// # Remarks
    /// The outer operation is always assumed to be adjointable, but does not
    /// need to be controllable in order for the combined operation to be
    /// controllable.
    ///
    /// # See Also
    /// - ApplyWithA
    /// - ApplyWithC
    /// - ApplyWithCA
    operation ApplyWith<'T>(outerOperation : ('T => Unit is Adj), innerOperation : ('T => Unit), target : 'T) : Unit {
        outerOperation(target);
        innerOperation(target);
        Adjoint outerOperation(target);
    }
    /// # Summary
    /// Given outer and inner operations, returns a new operation that
    /// conjugates the inner operation by the outer operation.
    ///
    /// # Input
    /// ## outerOperation
    /// The operation $U$ that should be used to conjugate $V$. Note that the
    /// outer operation $U$ needs to be adjointable, but does not
    /// need to be controllable.
    /// ## innerOperation
    /// The operation $V$ being conjugated.
    ///
    /// # Output
    /// A new operation whose action is represented by the unitary
    /// $U^{\dagger} V U$.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The type of the target on which each of the inner and outer operations
    /// act.
    ///
    /// # Remarks
    /// The outer operation is always assumed to be adjointable, but does not
    /// need to be controllable in order for the combined operation to be
    /// controllable.
    ///
    /// # See Also
    /// - ConjugatedByA
    /// - ConjugatedByC
    /// - ConjugatedByCA
    /// - ApplyWith
    function ConjugatedBy<'T>(outerOperation : ('T => Unit is Adj), innerOperation : ('T => Unit)) : ('T => Unit) {
        return ApplyWith(outerOperation, innerOperation, _);
    }

    /// # Summary
    /// Given two operations, applies one as conjugated with the other.
    ///
    /// # Description
    /// Given two operations, respectively described by unitary operators $U$
    /// and $V$, applies them in the sequence $U^{\dagger} V U$. That is,
    /// this operation implements the unitary operator given by $V$ conjugated
    /// with $U$.
    ///
    /// # Input
    /// ## outerOperation
    /// The operation $U$ that should be used to conjugate $V$. Note that the
    /// outer operation $U$ needs to be adjointable, but does not
    /// need to be controllable.
    /// ## innerOperation
    /// The operation $V$ being conjugated.
    /// ## target
    /// The input to be provided to the outer and inner operations.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The target on which each of the inner and outer operations act.
    ///
    /// # Remarks
    /// The outer operation is always assumed to be adjointable, but does not
    /// need to be controllable in order for the combined operation to be
    /// controllable.
    ///
    /// # See Also
    /// - ApplyWith
    /// - ApplyWithC
    /// - ApplyWithCA
    operation ApplyWithA<'T>(outerOperation : ('T => Unit is Adj), innerOperation : ('T => Unit is Adj), target : 'T) : Unit {
        body (...) {
            outerOperation(target);
            innerOperation(target);
            Adjoint outerOperation(target);
        }

        adjoint invert;
    }

    /// # Summary
    /// Given outer and inner operations, returns a new operation that
    /// conjugates the inner operation by the outer operation.
    ///
    /// # Input
    /// ## outerOperation
    /// The operation $U$ that should be used to conjugate $V$. Note that the
    /// outer operation $U$ needs to be adjointable, but does not
    /// need to be controllable.
    /// ## innerOperation
    /// The operation $V$ being conjugated.
    ///
    /// # Output
    /// A new operation whose action is represented by the unitary
    /// $U^{\dagger} V U$.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The type of the target on which each of the inner and outer operations
    /// act.
    ///
    /// # Remarks
    /// The outer operation is always assumed to be adjointable, but does not
    /// need to be controllable in order for the combined operation to be
    /// controllable.
    ///
    /// # See Also
    /// - ConjugatedByA
    /// - ConjugatedByC
    /// - ConjugatedByCA
    /// - ApplyWith
    function ConjugatedByA<'T>(outerOperation : ('T => Unit is Adj), innerOperation : ('T => Unit is Adj)) : ('T => Unit is Adj) {
        return ApplyWithA(outerOperation, innerOperation, _);
    }

    // # Summary
    /// Given two operations, applies one as conjugated with the other.
    ///
    /// # Description
    /// Given two operations, respectively described by unitary operators $U$
    /// and $V$, applies them in the sequence $U^{\dagger} V U$. That is,
    /// this operation implements the unitary operator given by $V$ conjugated
    /// with $U$.
    ///
    /// # Input
    /// ## outerOperation
    /// The operation $U$ that should be used to conjugate $V$. Note that the
    /// outer operation $U$ needs to be adjointable, but does not
    /// need to be controllable.
    /// ## innerOperation
    /// The operation $V$ being conjugated.
    /// ## target
    /// The input to be provided to the outer and inner operations.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The target on which each of the inner and outer operations act.
    ///
    /// # Remarks
    /// The outer operation is always assumed to be adjointable, but does not
    /// need to be controllable in order for the combined operation to be
    /// controllable.
    ///
    /// # See Also
    /// - ApplyWith
    /// - ApplyWithA
    /// - ApplyWithCA
    operation ApplyWithC<'T>(outerOperation : ('T => Unit is Adj), innerOperation : ('T => Unit is Ctl), target : 'T) : Unit {
        body (...) {
            outerOperation(target);
            innerOperation(target);
            Adjoint outerOperation(target);
        }

        controlled (controlRegister, ...) {
            outerOperation(target);
            Controlled innerOperation(controlRegister, target);
            Adjoint outerOperation(target);
        }
    }

    /// # Summary
    /// Given outer and inner operations, returns a new operation that
    /// conjugates the inner operation by the outer operation.
    ///
    /// # Input
    /// ## outerOperation
    /// The operation $U$ that should be used to conjugate $V$. Note that the
    /// outer operation $U$ needs to be adjointable, but does not
    /// need to be controllable.
    /// ## innerOperation
    /// The operation $V$ being conjugated.
    ///
    /// # Output
    /// A new operation whose action is represented by the unitary
    /// $U^{\dagger} V U$.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The type of the target on which each of the inner and outer operations
    /// act.
    ///
    /// # Remarks
    /// The outer operation is always assumed to be adjointable, but does not
    /// need to be controllable in order for the combined operation to be
    /// controllable.
    ///
    /// # See Also
    /// - ConjugatedBy
    /// - ConjugatedByA
    /// - ConjugatedByCA
    /// - ApplyWith
    function ConjugatedByC<'T>(outerOperation : ('T => Unit is Adj), innerOperation : ('T => Unit is Ctl)) : ('T => Unit is Ctl) {
        return ApplyWithC(outerOperation, innerOperation, _);
    }


    /// # Summary
    /// Given two operations, applies one as conjugated with the other.
    ///
    /// # Description
    /// Given two operations, respectively described by unitary operators $U$
    /// and $V$, applies them in the sequence $U^{\dagger} V U$. That is,
    /// this operation implements the unitary operator given by $V$ conjugated
    /// with $U$.
    ///
    /// # Input
    /// ## outerOperation
    /// The operation $U$ that should be used to conjugate $V$. Note that the
    /// outer operation $U$ needs to be adjointable, but does not
    /// need to be controllable.
    /// ## innerOperation
    /// The operation $V$ being conjugated.
    /// ## target
    /// The input to be provided to the outer and inner operations.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The target on which each of the inner and outer operations act.
    ///
    /// # Remarks
    /// The outer operation is always assumed to be adjointable, but does not
    /// need to be controllable in order for the combined operation to be
    /// controllable.
    ///
    /// # See Also
    /// - ApplyWith
    /// - ApplyWithA
    /// - ApplyWithC
    operation ApplyWithCA<'T>(outerOperation : ('T => Unit is Adj), innerOperation : ('T => Unit is Adj + Ctl), target : 'T) : Unit {
        body (...) {
            outerOperation(target);
            innerOperation(target);
            Adjoint outerOperation(target);
        }

        adjoint auto;

        controlled (controlRegister, ...) {
            outerOperation(target);
            Controlled innerOperation(controlRegister, target);
            Adjoint outerOperation(target);
        }

        controlled adjoint auto;
    }


    /// # Summary
    /// Given outer and inner operations, returns a new operation that
    /// conjugates the inner operation by the outer operation.
    ///
    /// # Input
    /// ## outerOperation
    /// The operation $U$ that should be used to conjugate $V$. Note that the
    /// outer operation $U$ needs to be adjointable, but does not
    /// need to be controllable.
    /// ## innerOperation
    /// The operation $V$ being conjugated.
    ///
    /// # Output
    /// A new operation whose action is represented by the unitary
    /// $U^{\dagger} V U$.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The type of the target on which each of the inner and outer operations
    /// act.
    ///
    /// # Remarks
    /// The outer operation is always assumed to be adjointable, but does not
    /// need to be controllable in order for the combined operation to be
    /// controllable.
    ///
    /// # See Also
    /// - ConjugatedByA
    /// - ConjugatedByC
    /// - ConjugatedByCA
    /// - ApplyWith
    function ConjugatedByCA<'T>(outerOperation : ('T => Unit is Adj), innerOperation : ('T => Unit is Adj + Ctl)) : ('T => Unit is Adj + Ctl) {
        return ApplyWithCA(outerOperation, innerOperation, _);
    }

}
