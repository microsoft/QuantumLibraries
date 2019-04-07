// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {

    // #region No functors

    /// # Summary
    /// Given an operation that accepts some input, a function that
    /// returns an output compatible with that operation, and an input to that
    /// function, applies the operation using the function to transform the
    /// input to a form expected by the operation.
    ///
    /// # Input
    /// ## fn
    /// A function that transforms the given input into a form expected by the
    /// operation.
    /// ## op
    /// The operation to be applied.
    /// ## input
    /// An input to the function.
    ///
    /// # See Also
    /// - ApplyWithInputTransformationA
    /// - ApplyWithInputTransformationC
    /// - ApplyWithInputTransformationCA
    /// - TransformedOperation
    ///
    /// # Example
    /// The following call uses
    /// @"Microsoft.Quantum.Arithmetic.LittleEndianAsBigEndian" to apply
    /// an operation designed for
    /// @"Microsoft.Quantum.Arithmetic.BigEndian" inputs to an input of type
    /// @"Microsoft.Quantum.Arithmetic.LittleEndian":
    /// ```qsharp
    /// ApplyWithInputTransformation(LittleEndianAsBigEndian, ApplyXorInPlaceBE, LittleEndian(qubits));
    /// ```
    operation ApplyWithInputTransformation<'T, 'U>(fn : ('U -> 'T), op : ('T => Unit), input : 'U) : Unit {
        op(fn(input));
    }

    /// # Summary
    /// Given a function and an operation, returns a new operation whose
    /// input is transformed by the given function.
    ///
    /// # Input
    /// ## fn
    /// A function that transforms the given input into a form expected by the
    /// operation.
    /// ## op
    /// The operation to be transformed.
    ///
    /// # Output
    /// A new operation tbat calls `fn` with its input, then passes the
    /// resulting output to `op`.
    ///
    /// # See Also
    /// - TransformedOperationA
    /// - TransformedOperationC
    /// - TransformedOperationCA
    /// - ApplyWithInputTransformation
    /// - Compose
    ///
    /// # Example
    /// The following call uses
    /// @"Microsoft.Quantum.Arithmetic.LittleEndianAsBigEndian" to transform
    /// an operation designed for
    /// @"Microsoft.Quantum.Arithmetic.BigEndian" inputs into an operation
    /// that accepts inputs of type
    /// @"Microsoft.Quantum.Arithmetic.LittleEndian":
    /// ```qsharp
    /// let leOp = TransformedOperation(LittleEndianAsBigEndian, ApplyXorInPlaceBE);
    /// ```
    function TransformedOperation<'T, 'U>(fn : ('U -> 'T), op : ('T => Unit)) : ('U => Unit) {
        return ApplyWithInputTransformation(fn, op, _);
    }

    // #endregion


    // #region Adjoint

    /// # Summary
    /// Given an operation that accepts some input, a function that
    /// returns an output compatible with that operation, and an input to that
    /// function, applies the operation using the function to transform the
    /// input to a form expected by the operation.
    ///
    /// # Input
    /// ## fn
    /// A function that transforms the given input into a form expected by the
    /// operation.
    /// ## op
    /// The operation to be applied.
    /// ## input
    /// An input to the function.
    ///
    /// # See Also
    /// - ApplyWithInputTransformation
    /// - ApplyWithInputTransformationC
    /// - ApplyWithInputTransformationCA
    /// - TransformedOperation
    ///
    /// # Example
    /// The following call uses
    /// @"Microsoft.Quantum.Arithmetic.LittleEndianAsBigEndian" to apply
    /// an operation designed for
    /// @"Microsoft.Quantum.Arithmetic.BigEndian" inputs to an input of type
    /// @"Microsoft.Quantum.Arithmetic.LittleEndian":
    /// ```qsharp
    /// ApplyWithInputTransformation(LittleEndianAsBigEndian, ApplyXorInPlaceBE, LittleEndian(qubits));
    /// ```
    operation ApplyWithInputTransformationA<'T, 'U>(fn : ('U -> 'T), op : ('T => Unit : Adjoint), input : 'U) : Unit {
        body (...) { op(fn(input)); }
        adjoint auto;
    }

    /// # Summary
    /// Given a function and an operation, returns a new operation whose
    /// input is transformed by the given function.
    ///
    /// # Input
    /// ## fn
    /// A function that transforms the given input into a form expected by the
    /// operation.
    /// ## op
    /// The operation to be transformed.
    ///
    /// # Output
    /// A new operation tbat calls `fn` with its input, then passes the
    /// resulting output to `op`.
    ///
    /// # See Also
    /// - TransformedOperation
    /// - TransformedOperationC
    /// - TransformedOperationCA
    /// - ApplyWithInputTransformation
    /// - Compose
    ///
    /// # Example
    /// The following call uses
    /// @"Microsoft.Quantum.Arithmetic.LittleEndianAsBigEndian" to transform
    /// an operation designed for
    /// @"Microsoft.Quantum.Arithmetic.BigEndian" inputs into an operation
    /// that accepts inputs of type
    /// @"Microsoft.Quantum.Arithmetic.LittleEndian":
    /// ```qsharp
    /// let leOp = TransformedOperation(LittleEndianAsBigEndian, ApplyXorInPlaceBE);
    /// ```
    function TransformedOperationA<'T, 'U>(fn : ('U -> 'T), op : ('T => Unit : Adjoint)) : ('U => Unit : Adjoint) {
        return ApplyWithInputTransformationA(fn, op, _);
    }

    // #endregion



    // #region Adjoint

    /// # Summary
    /// Given an operation that accepts some input, a function that
    /// returns an output compatible with that operation, and an input to that
    /// function, applies the operation using the function to transform the
    /// input to a form expected by the operation.
    ///
    /// # Input
    /// ## fn
    /// A function that transforms the given input into a form expected by the
    /// operation.
    /// ## op
    /// The operation to be applied.
    /// ## input
    /// An input to the function.
    ///
    /// # See Also
    /// - ApplyWithInputTransformation
    /// - ApplyWithInputTransformationA
    /// - ApplyWithInputTransformationCA
    /// - TransformedOperation
    ///
    /// # Example
    /// The following call uses
    /// @"Microsoft.Quantum.Arithmetic.LittleEndianAsBigEndian" to apply
    /// an operation designed for
    /// @"Microsoft.Quantum.Arithmetic.BigEndian" inputs to an input of type
    /// @"Microsoft.Quantum.Arithmetic.LittleEndian":
    /// ```qsharp
    /// ApplyWithInputTransformation(LittleEndianAsBigEndian, ApplyXorInPlaceBE, LittleEndian(qubits));
    /// ```
    operation ApplyWithInputTransformationC<'T, 'U>(fn : ('U -> 'T), op : ('T => Unit : Controlled), input : 'U) : Unit {
        body (...) { op(fn(input)); }
        controlled auto;
    }

    /// # Summary
    /// Given a function and an operation, returns a new operation whose
    /// input is transformed by the given function.
    ///
    /// # Input
    /// ## fn
    /// A function that transforms the given input into a form expected by the
    /// operation.
    /// ## op
    /// The operation to be transformed.
    ///
    /// # Output
    /// A new operation tbat calls `fn` with its input, then passes the
    /// resulting output to `op`.
    ///
    /// # See Also
    /// - TransformedOperation
    /// - TransformedOperationA
    /// - TransformedOperationCA
    /// - ApplyWithInputTransformation
    /// - Compose
    ///
    /// # Example
    /// The following call uses
    /// @"Microsoft.Quantum.Arithmetic.LittleEndianAsBigEndian" to transform
    /// an operation designed for
    /// @"Microsoft.Quantum.Arithmetic.BigEndian" inputs into an operation
    /// that accepts inputs of type
    /// @"Microsoft.Quantum.Arithmetic.LittleEndian":
    /// ```qsharp
    /// let leOp = TransformedOperation(LittleEndianAsBigEndian, ApplyXorInPlaceBE);
    /// ```
    function TransformedOperationC<'T, 'U>(fn : ('U -> 'T), op : ('T => Unit : Controlled)) : ('U => Unit : Controlled) {
        return ApplyWithInputTransformationC(fn, op, _);
    }

    // #endregion

    
    // #region Adjoint

    /// # Summary
    /// Given an operation that accepts some input, a function that
    /// returns an output compatible with that operation, and an input to that
    /// function, applies the operation using the function to transform the
    /// input to a form expected by the operation.
    ///
    /// # Input
    /// ## fn
    /// A function that transforms the given input into a form expected by the
    /// operation.
    /// ## op
    /// The operation to be applied.
    /// ## input
    /// An input to the function.
    ///
    /// # See Also
    /// - ApplyWithInputTransformation
    /// - ApplyWithInputTransformationA
    /// - ApplyWithInputTransformationC
    /// - TransformedOperation
    ///
    /// # Example
    /// The following call uses
    /// @"Microsoft.Quantum.Arithmetic.LittleEndianAsBigEndian" to apply
    /// an operation designed for
    /// @"Microsoft.Quantum.Arithmetic.BigEndian" inputs to an input of type
    /// @"Microsoft.Quantum.Arithmetic.LittleEndian":
    /// ```qsharp
    /// ApplyWithInputTransformation(LittleEndianAsBigEndian, ApplyXorInPlaceBE, LittleEndian(qubits));
    /// ```
    operation ApplyWithInputTransformationCA<'T, 'U>(fn : ('U -> 'T), op : ('T => Unit : Adjoint, Controlled), input : 'U) : Unit {
        body (...) { op(fn(input)); }
        adjoint auto;
        controlled auto;
        adjoint controlled auto;
    }

    /// # Summary
    /// Given a function and an operation, returns a new operation whose
    /// input is transformed by the given function.
    ///
    /// # Input
    /// ## fn
    /// A function that transforms the given input into a form expected by the
    /// operation.
    /// ## op
    /// The operation to be transformed.
    ///
    /// # Output
    /// A new operation tbat calls `fn` with its input, then passes the
    /// resulting output to `op`.
    ///
    /// # See Also
    /// - TransformedOperation
    /// - TransformedOperationA
    /// - TransformedOperationCA
    /// - ApplyWithInputTransformation
    /// - Compose
    ///
    /// # Example
    /// The following call uses
    /// @"Microsoft.Quantum.Arithmetic.LittleEndianAsBigEndian" to transform
    /// an operation designed for
    /// @"Microsoft.Quantum.Arithmetic.BigEndian" inputs into an operation
    /// that accepts inputs of type
    /// @"Microsoft.Quantum.Arithmetic.LittleEndian":
    /// ```qsharp
    /// let leOp = TransformedOperation(LittleEndianAsBigEndian, ApplyXorInPlaceBE);
    /// ```
    function TransformedOperationCA<'T, 'U>(fn : ('U -> 'T), op : ('T => Unit : Controlled, Adjoint)) : ('U => Unit : Controlled, Adjoint) {
        return ApplyWithInputTransformationCA(fn, op, _);
    }

    // #endregion

}
