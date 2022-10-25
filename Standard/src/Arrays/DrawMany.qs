// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Arrays {

    /// # Summary
    /// Repeats an operation for a given number of samples, collecting its outputs
    /// in an array.
    ///
    /// # Input
    /// ## op
    /// The operation to be called repeatedly.
    /// ## nSamples
    /// The number of samples of calling `op` to collect.
    /// ## input
    /// The input to be passed to `op`.
    ///
    /// # Type Parameters
    /// ## TInput
    /// The type of input expected by `op`.
    /// ## TOutput
    /// The type of output returned by `op`.
    ///
    /// # Example
    /// The following samples an integer, given an operation
    /// that samples a single bit at a time.
    /// ```qsharp
    /// let randomInteger = BoolArrayAsInt(DrawMany(SampleRandomBit, 16, ()));
    /// ```
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.Repeat
    operation DrawMany<'TInput, 'TOutput>(op : ('TInput => 'TOutput), nSamples : Int, input : 'TInput)
    : 'TOutput[] {
        if nSamples == 0 {
            return [];
        }

        let first = op(input);
        mutable outputs = [first, size = nSamples];
        for idx in 1..nSamples - 1 {
            set outputs w/= idx <- op(input);
        }
        return outputs;
    }

}
