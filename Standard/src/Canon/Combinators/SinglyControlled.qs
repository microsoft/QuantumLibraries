// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {

    /// # Summary
    /// Given a controllable operation, returns a controlled version of that operation
    /// accepting exactly one control qubit.
    ///
    /// # Input
    /// ## op
    /// The operation to be controlled.
    ///
    /// # Output
    /// A controlled variant of `op` accepting exactly one control qubit.
    ///
    /// # Example
    /// To add the weight (number of "1" bits) of a control register to
    /// a target register:
    /// ```qsharp
    /// ApplyToEachCA(
    ///     SinglyControlled(IncrementByInteger)(_, (1, target)),
    ///     controls)
    /// );
    /// ```
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.SinglyControlledA
    function SinglyControlled<'T>(op : 'T => Unit is Ctl) : (Qubit, 'T) => Unit is Ctl {
        return (ctrl, originalInput) => Controlled op([ctrl], originalInput);
    }

    /// # Summary
    /// Given a controllable operation, returns a controlled version of that operation
    /// accepting exactly one control qubit.
    ///
    /// # Input
    /// ## op
    /// The operation to be controlled.
    ///
    /// # Output
    /// A controlled variant of `op` accepting exactly one control qubit.
    ///
    /// # Example
    /// To add the weight (number of "1" bits) of a control register to
    /// a target register:
    /// ```qsharp
    /// ApplyToEachCA(
    ///     SinglyControlledA(IncrementByInteger)(_, (1, target)),
    ///     controls)
    /// );
    /// ```
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.SinglyControlled
    function SinglyControlledA<'T>(op : 'T => Unit is Adj + Ctl) : (Qubit, 'T) => Unit is Adj + Ctl {
        return (ctrl, originalInput) => Controlled op([ctrl], originalInput);
    }

}
