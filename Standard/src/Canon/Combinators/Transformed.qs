// Copyright (c) Microsoft Corporation.
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
    /// - Microsoft.Quantum.Canon.ApplyWithInputTransformationA
    /// - Microsoft.Quantum.Canon.ApplyWithInputTransformationC
    /// - Microsoft.Quantum.Canon.ApplyWithInputTransformationCA
    /// - Microsoft.Quantum.Canon.TransformedOperation
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
    /// A new operation that calls `fn` with its input, then passes the
    /// resulting output to `op`.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.TransformedOperationA
    /// - Microsoft.Quantum.Canon.TransformedOperationC
    /// - Microsoft.Quantum.Canon.TransformedOperationCA
    /// - Microsoft.Quantum.Canon.ApplyWithInputTransformation
    /// - Microsoft.Quantum.Canon.Composed
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
    /// - Microsoft.Quantum.Canon.ApplyWithInputTransformation
    /// - Microsoft.Quantum.Canon.ApplyWithInputTransformationC
    /// - Microsoft.Quantum.Canon.ApplyWithInputTransformationCA
    /// - Microsoft.Quantum.Canon.TransformedOperation
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
    operation ApplyWithInputTransformationA<'T, 'U>(fn : ('U -> 'T), op : ('T => Unit is Adj), input : 'U) : Unit {
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
    /// A new operation that calls `fn` with its input, then passes the
    /// resulting output to `op`.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.TransformedOperation
    /// - Microsoft.Quantum.Canon.TransformedOperationC
    /// - Microsoft.Quantum.Canon.TransformedOperationCA
    /// - Microsoft.Quantum.Canon.ApplyWithInputTransformation
    /// - Microsoft.Quantum.Canon.Composed
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
    function TransformedOperationA<'T, 'U>(fn : ('U -> 'T), op : ('T => Unit is Adj)) : ('U => Unit is Adj) {
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
    /// - Microsoft.Quantum.Canon.ApplyWithInputTransformation
    /// - Microsoft.Quantum.Canon.ApplyWithInputTransformationA
    /// - Microsoft.Quantum.Canon.ApplyWithInputTransformationCA
    /// - Microsoft.Quantum.Canon.TransformedOperation
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
    operation ApplyWithInputTransformationC<'T, 'U>(fn : ('U -> 'T), op : ('T => Unit is Ctl), input : 'U) : Unit {
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
    /// A new operation that calls `fn` with its input, then passes the
    /// resulting output to `op`.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.TransformedOperation
    /// - Microsoft.Quantum.Canon.TransformedOperationA
    /// - Microsoft.Quantum.Canon.TransformedOperationCA
    /// - Microsoft.Quantum.Canon.ApplyWithInputTransformation
    /// - Microsoft.Quantum.Canon.Composed
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
    function TransformedOperationC<'T, 'U>(fn : ('U -> 'T), op : ('T => Unit is Ctl)) : ('U => Unit is Ctl) {
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
    /// - Microsoft.Quantum.Canon.ApplyWithInputTransformation
    /// - Microsoft.Quantum.Canon.ApplyWithInputTransformationA
    /// - Microsoft.Quantum.Canon.ApplyWithInputTransformationC
    /// - Microsoft.Quantum.Canon.TransformedOperation
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
    operation ApplyWithInputTransformationCA<'T, 'U>(fn : ('U -> 'T), op : ('T => Unit is Adj + Ctl), input : 'U) : Unit {
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
    /// A new operation that calls `fn` with its input, then passes the
    /// resulting output to `op`.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.TransformedOperation
    /// - Microsoft.Quantum.Canon.TransformedOperationA
    /// - Microsoft.Quantum.Canon.TransformedOperationCA
    /// - Microsoft.Quantum.Canon.ApplyWithInputTransformation
    /// - Microsoft.Quantum.Canon.Composed
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
    function TransformedOperationCA<'T, 'U>(fn : ('U -> 'T), op : ('T => Unit is Adj + Ctl)) : ('U => Unit is Adj + Ctl) {
        return ApplyWithInputTransformationCA(fn, op, _);
    }

    // #endregion

}
