// Copyright (c) Microsoft Corporation.
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
    /// - Microsoft.Quantum.Canon.ApplyWithA
    /// - Microsoft.Quantum.Canon.ApplyWithC
    /// - Microsoft.Quantum.Canon.ApplyWithCA
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
    /// - Microsoft.Quantum.Canon.ConjugatedByA
    /// - Microsoft.Quantum.Canon.ConjugatedByC
    /// - Microsoft.Quantum.Canon.ConjugatedByCA
    /// - Microsoft.Quantum.Canon.ApplyWith
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
    /// - Microsoft.Quantum.Canon.ApplyWith
    /// - Microsoft.Quantum.Canon.ApplyWithC
    /// - Microsoft.Quantum.Canon.ApplyWithCA
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
    /// - Microsoft.Quantum.Canon.ConjugatedByA
    /// - Microsoft.Quantum.Canon.ConjugatedByC
    /// - Microsoft.Quantum.Canon.ConjugatedByCA
    /// - Microsoft.Quantum.Canon.ApplyWith
    function ConjugatedByA<'T>(outerOperation : ('T => Unit is Adj), innerOperation : ('T => Unit is Adj)) : ('T => Unit is Adj) {
        return ApplyWithA(outerOperation, innerOperation, _);
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
    /// - Microsoft.Quantum.Canon.ApplyWith
    /// - Microsoft.Quantum.Canon.ApplyWithA
    /// - Microsoft.Quantum.Canon.ApplyWithCA
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
    /// - Microsoft.Quantum.Canon.ConjugatedBy
    /// - Microsoft.Quantum.Canon.ConjugatedByA
    /// - Microsoft.Quantum.Canon.ConjugatedByCA
    /// - Microsoft.Quantum.Canon.ApplyWith
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
    /// - Microsoft.Quantum.Canon.ApplyWith
    /// - Microsoft.Quantum.Canon.ApplyWithA
    /// - Microsoft.Quantum.Canon.ApplyWithC
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
    /// - Microsoft.Quantum.Canon.ConjugatedByA
    /// - Microsoft.Quantum.Canon.ConjugatedByC
    /// - Microsoft.Quantum.Canon.ConjugatedByCA
    /// - Microsoft.Quantum.Canon.ApplyWith
    function ConjugatedByCA<'T>(outerOperation : ('T => Unit is Adj), innerOperation : ('T => Unit is Adj + Ctl)) : ('T => Unit is Adj + Ctl) {
        return ApplyWithCA(outerOperation, innerOperation, _);
    }

}
