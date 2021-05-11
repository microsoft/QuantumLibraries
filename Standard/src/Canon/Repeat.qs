// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {

    /// # Summary
    /// Repeats an operation a given number of times.
    ///
    /// # Input
    /// ## op
    /// The operation to be called repeatedly.
    /// ## nTimes
    /// The number of times to call `op`.
    /// ## input
    /// The input to be passed to `op`.
    ///
    /// # Type Parameters
    /// ## TInput
    /// The type of input expected by `op`.
    ///
    /// # Example
    /// The following are equivalent:
    /// ```qsharp
    /// Repeat(U, 17, target);
    /// (Bound(ConstantArray(17, U)))(target);
    /// ```
    ///
    /// # See Also
    /// - Microsoft.Quantum.Arrays.DrawMany
    /// - Microsoft.Quantum.Canon.RepeatA
    /// - Microsoft.Quantum.Canon.RepeatC
    /// - Microsoft.Quantum.Canon.RepeatCA
    operation Repeat<'TInput>(op : ('TInput => Unit), nTimes : Int, input : 'TInput) : Unit {
        for idx in 0..nTimes - 1 {
            op(input);
        }
    }

    /// # Summary
    /// Repeats an operation a given number of times.
    ///
    /// # Input
    /// ## op
    /// The operation to be called repeatedly.
    /// ## nTimes
    /// The number of times to call `op`.
    /// ## input
    /// The input to be passed to `op`.
    ///
    /// # Type Parameters
    /// ## TInput
    /// The type of input expected by `op`.
    ///
    /// # Example
    /// The following are equivalent:
    /// ```qsharp
    /// RepeatA(U, 17, target);
    /// (BoundA(ConstantArray(17, U)))(target);
    /// ```
    ///
    /// # See Also
    /// - Microsoft.Quantum.Arrays.DrawMany
    /// - Microsoft.Quantum.Canon.Repeat
    /// - Microsoft.Quantum.Canon.RepeatC
    /// - Microsoft.Quantum.Canon.RepeatCA
    operation RepeatA<'TInput>(op : ('TInput => Unit is Adj), nTimes : Int, input : 'TInput) : Unit is Adj {
        for idx in 0..nTimes - 1 {
            op(input);
        }
    }
    
    /// # Summary
    /// Repeats an operation a given number of times.
    ///
    /// # Input
    /// ## op
    /// The operation to be called repeatedly.
    /// ## nTimes
    /// The number of times to call `op`.
    /// ## input
    /// The input to be passed to `op`.
    ///
    /// # Type Parameters
    /// ## TInput
    /// The type of input expected by `op`.
    ///
    /// # Example
    /// The following are equivalent:
    /// ```qsharp
    /// RepeatC(U, 17, target);
    /// (BoundC(ConstantArray(17, U)))(target);
    /// ```
    ///
    /// # See Also
    /// - Microsoft.Quantum.Arrays.DrawMany
    /// - Microsoft.Quantum.Canon.Repeat
    /// - Microsoft.Quantum.Canon.RepeatA
    /// - Microsoft.Quantum.Canon.RepeatCA
    operation RepeatC<'TInput>(op : ('TInput => Unit is Ctl), nTimes : Int, input : 'TInput) : Unit is Ctl {
        for idx in 0..nTimes - 1 {
            op(input);
        }
    }

    /// # Summary
    /// Repeats an operation a given number of times.
    ///
    /// # Input
    /// ## op
    /// The operation to be called repeatedly.
    /// ## nTimes
    /// The number of times to call `op`.
    /// ## input
    /// The input to be passed to `op`.
    ///
    /// # Type Parameters
    /// ## TInput
    /// The type of input expected by `op`.
    ///
    /// # Example
    /// The following are equivalent:
    /// ```qsharp
    /// RepeatCA(U, 17, target);
    /// (BoundCA(ConstantArray(17, U)))(target);
    /// ```
    ///
    /// # See Also
    /// - Microsoft.Quantum.Arrays.DrawMany
    /// - Microsoft.Quantum.Canon.Repeat
    /// - Microsoft.Quantum.Canon.RepeatA
    /// - Microsoft.Quantum.Canon.RepeatC
    operation RepeatCA<'TInput>(op : ('TInput => Unit is Adj + Ctl), nTimes : Int, input : 'TInput) : Unit is Adj + Ctl { 
        for idx in 0..nTimes - 1 {
            op(input);
        }
    }

}
