// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {
    open Microsoft.Quantum.Arrays;

    /// # Summary
    /// Applies an operation to the first element of an array.
    ///
    /// # Description
    /// Given an operation `op` and an array of targets `targets`,
    /// applies `op(Head(targets))`.
    ///
    /// # Input
    /// ## op
    /// An operation to be applied.
    /// ## targets
    /// An array of targets, of which the first will be applied to `op`.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The input type of the operation to be applied.
    ///
    /// # Example
    /// The following Q# snippets are equivalent:
    /// ```qsharp
    /// ApplyToHead(H, register);
    /// H(Head(register));
    /// ```
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApplyToHeadA
    /// - Microsoft.Quantum.Canon.ApplyToHeadC
    /// - Microsoft.Quantum.Canon.ApplyToHeadCA
    operation ApplyToHead<'T>(op : ('T => Unit), targets : 'T[]) : Unit {
        op(Head(targets));
    }

    /// # Summary
    /// Applies an operation to the first element of an array.
    ///
    /// # Description
    /// Given an operation `op` and an array of targets `targets`,
    /// applies `op(Head(targets))`.
    ///
    /// # Input
    /// ## op
    /// An operation to be applied.
    /// ## targets
    /// An array of targets, of which the first will be applied to `op`.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The input type of the operation to be applied.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApplyToHead
    /// - Microsoft.Quantum.Canon.ApplyToHeadC
    /// - Microsoft.Quantum.Canon.ApplyToHeadCA
    operation ApplyToHeadA<'T>(op : ('T => Unit is Adj), targets : 'T[]) : Unit is Adj {
        op(Head(targets));
    }

    /// # Summary
    /// Applies an operation to the first element of an array.
    ///
    /// # Description
    /// Given an operation `op` and an array of targets `targets`,
    /// applies `op(Head(targets))`.
    ///
    /// # Input
    /// ## op
    /// An operation to be applied.
    /// ## targets
    /// An array of targets, of which the first will be applied to `op`.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The input type of the operation to be applied.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApplyToHead
    /// - Microsoft.Quantum.Canon.ApplyToHeadA
    /// - Microsoft.Quantum.Canon.ApplyToHeadCA
    operation ApplyToHeadC<'T>(op : ('T => Unit is Ctl), targets : 'T[]) : Unit is Ctl {
        op(Head(targets));
    }

    /// # Summary
    /// Applies an operation to the first element of an array.
    ///
    /// # Description
    /// Given an operation `op` and an array of targets `targets`,
    /// applies `op(Head(targets))`.
    ///
    /// # Input
    /// ## op
    /// An operation to be applied.
    /// ## targets
    /// An array of targets, of which the first will be applied to `op`.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The input type of the operation to be applied.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApplyToHead
    /// - Microsoft.Quantum.Canon.ApplyToHeadA
    /// - Microsoft.Quantum.Canon.ApplyToHeadC
    operation ApplyToHeadCA<'T>(op : ('T => Unit is Adj+Ctl), targets : 'T[]) : Unit is Adj+Ctl {
        op(Head(targets));
    }

    /// # Summary
    /// Applies an operation to all but the first element of an array.
    ///
    /// # Description
    /// Given an operation `op` and an array of targets `targets`,
    /// applies `op(Rest(targets))`.
    ///
    /// # Input
    /// ## op
    /// An operation to be applied.
    /// ## targets
    /// An array of targets, of which all but the first will be applied to `op`.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The input type of the operation to be applied.
    ///
    /// # Example
    /// The following Q# snippets are equivalent:
    /// ```qsharp
    /// ApplyToRest(ApplyCNOTChain, register);
    /// ApplyCNOTChain(Rest(register));
    /// ```
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApplyToRestA
    /// - Microsoft.Quantum.Canon.ApplyToRestC
    /// - Microsoft.Quantum.Canon.ApplyToRestCA
    operation ApplyToRest<'T>(op : ('T[] => Unit), targets : 'T[]) : Unit {
        op(Rest(targets));
    }

    /// # Summary
    /// Applies an operation to all but the first element of an array.
    ///
    /// # Description
    /// Given an operation `op` and an array of targets `targets`,
    /// applies `op(Rest(targets))`.
    ///
    /// # Input
    /// ## op
    /// An operation to be applied.
    /// ## targets
    /// An array of targets, of which all but the first will be applied to `op`.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The input type of the operation to be applied.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApplyToRest
    /// - Microsoft.Quantum.Canon.ApplyToRestC
    /// - Microsoft.Quantum.Canon.ApplyToRestCA
    operation ApplyToRestA<'T>(op : ('T[] => Unit is Adj), targets : 'T[]) : Unit is Adj {
        op(Rest(targets));
    }

    /// # Summary
    /// Applies an operation to all but the first element of an array.
    ///
    /// # Description
    /// Given an operation `op` and an array of targets `targets`,
    /// applies `op(Rest(targets))`.
    ///
    /// # Input
    /// ## op
    /// An operation to be applied.
    /// ## targets
    /// An array of targets, of which all but the first will be applied to `op`.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The input type of the operation to be applied.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApplyToRest
    /// - Microsoft.Quantum.Canon.ApplyToRestA
    /// - Microsoft.Quantum.Canon.ApplyToRestCA
    operation ApplyToRestC<'T>(op : ('T[] => Unit is Ctl), targets : 'T[]) : Unit is Ctl {
        op(Rest(targets));
    }

    /// # Summary
    /// Applies an operation to all but the first element of an array.
    ///
    /// # Description
    /// Given an operation `op` and an array of targets `targets`,
    /// applies `op(Rest(targets))`.
    ///
    /// # Input
    /// ## op
    /// An operation to be applied.
    /// ## targets
    /// An array of targets, of which all but the first will be applied to `op`.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The input type of the operation to be applied.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApplyToRest
    /// - Microsoft.Quantum.Canon.ApplyToRestA
    /// - Microsoft.Quantum.Canon.ApplyToRestC
    operation ApplyToRestCA<'T>(op : ('T[] => Unit is Adj+Ctl), targets : 'T[]) : Unit is Adj+Ctl {
        op(Rest(targets));
    }

    /// # Summary
    /// Applies an operation to the last element of an array.
    ///
    /// # Description
    /// Given an operation `op` and an array of targets `targets`,
    /// applies `op(Tail(targets))`.
    ///
    /// # Input
    /// ## op
    /// An operation to be applied.
    /// ## targets
    /// An array of targets, of which the last will be applied to `op`.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The input type of the operation to be applied.
    ///
    /// # Example
    /// The following Q# snippets are equivalent:
    /// ```qsharp
    /// ApplyToTail(H, register);
    /// H(Tail(register));
    /// ```
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApplyToTailA
    /// - Microsoft.Quantum.Canon.ApplyToTailC
    /// - Microsoft.Quantum.Canon.ApplyToTailCA
    operation ApplyToTail<'T>(op : ('T => Unit), targets : 'T[]) : Unit {
        op(Tail(targets));
    }

    /// # Summary
    /// Applies an operation to the last element of an array.
    ///
    /// # Description
    /// Given an operation `op` and an array of targets `targets`,
    /// applies `op(Tail(targets))`.
    ///
    /// # Input
    /// ## op
    /// An operation to be applied.
    /// ## targets
    /// An array of targets, of which the last will be applied to `op`.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The input type of the operation to be applied.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApplyToTail
    /// - Microsoft.Quantum.Canon.ApplyToTailC
    /// - Microsoft.Quantum.Canon.ApplyToTailCA
    operation ApplyToTailA<'T>(op : ('T => Unit is Adj), targets : 'T[]) : Unit is Adj {
        op(Tail(targets));
    }

    /// # Summary
    /// Applies an operation to the last element of an array.
    ///
    /// # Description
    /// Given an operation `op` and an array of targets `targets`,
    /// applies `op(Tail(targets))`.
    ///
    /// # Input
    /// ## op
    /// An operation to be applied.
    /// ## targets
    /// An array of targets, of which the last will be applied to `op`.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The input type of the operation to be applied.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApplyToTail
    /// - Microsoft.Quantum.Canon.ApplyToTailA
    /// - Microsoft.Quantum.Canon.ApplyToTailCA
    operation ApplyToTailC<'T>(op : ('T => Unit is Ctl), targets : 'T[]) : Unit is Ctl {
        op(Tail(targets));
    }

    /// # Summary
    /// Applies an operation to the last element of an array.
    ///
    /// # Description
    /// Given an operation `op` and an array of targets `targets`,
    /// applies `op(Tail(targets))`.
    ///
    /// # Input
    /// ## op
    /// An operation to be applied.
    /// ## targets
    /// An array of targets, of which the last will be applied to `op`.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The input type of the operation to be applied.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApplyToTail
    /// - Microsoft.Quantum.Canon.ApplyToTailA
    /// - Microsoft.Quantum.Canon.ApplyToTailC
    operation ApplyToTailCA<'T>(op : ('T => Unit is Adj+Ctl), targets : 'T[]) : Unit is Adj+Ctl {
        op(Tail(targets));
    }

    /// # Summary
    /// Applies an operation to all but the last element of an array.
    ///
    /// # Description
    /// Given an operation `op` and an array of targets `targets`,
    /// applies `op(Most(targets))`.
    ///
    /// # Input
    /// ## op
    /// An operation to be applied.
    /// ## targets
    /// An array of targets, of which all but the last will be applied to `op`.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The input type of the operation to be applied.
    ///
    /// # Example
    /// The following Q# snippets are equivalent:
    /// ```qsharp
    /// ApplyToMost(ApplyCNOTChain, register);
    /// ApplyCNOTChain(Most(register));
    /// ```
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApplyToMostA
    /// - Microsoft.Quantum.Canon.ApplyToMostC
    /// - Microsoft.Quantum.Canon.ApplyToMostCA
    operation ApplyToMost<'T>(op : ('T[] => Unit), targets : 'T[]) : Unit {
        op(Most(targets));
    }

    /// # Summary
    /// Applies an operation to all but the last element of an array.
    ///
    /// # Description
    /// Given an operation `op` and an array of targets `targets`,
    /// applies `op(Most(targets))`.
    ///
    /// # Input
    /// ## op
    /// An operation to be applied.
    /// ## targets
    /// An array of targets, of which all but the last will be applied to `op`.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The input type of the operation to be applied.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApplyToMost
    /// - Microsoft.Quantum.Canon.ApplyToMostC
    /// - Microsoft.Quantum.Canon.ApplyToMostCA
    operation ApplyToMostA<'T>(op : ('T[] => Unit is Adj), targets : 'T[]) : Unit is Adj {
        op(Most(targets));
    }

    /// # Summary
    /// Applies an operation to all but the last element of an array.
    ///
    /// # Description
    /// Given an operation `op` and an array of targets `targets`,
    /// applies `op(Most(targets))`.
    ///
    /// # Input
    /// ## op
    /// An operation to be applied.
    /// ## targets
    /// An array of targets, of which all but the last will be applied to `op`.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The input type of the operation to be applied.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApplyToMost
    /// - Microsoft.Quantum.Canon.ApplyToMostA
    /// - Microsoft.Quantum.Canon.ApplyToMostCA
    operation ApplyToMostC<'T>(op : ('T[] => Unit is Ctl), targets : 'T[]) : Unit is Ctl {
        op(Most(targets));
    }

    /// # Summary
    /// Applies an operation to all but the last element of an array.
    ///
    /// # Description
    /// Given an operation `op` and an array of targets `targets`,
    /// applies `op(Most(targets))`.
    ///
    /// # Input
    /// ## op
    /// An operation to be applied.
    /// ## targets
    /// An array of targets, of which all but the last will be applied to `op`.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The input type of the operation to be applied.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApplyToMost
    /// - Microsoft.Quantum.Canon.ApplyToMostA
    /// - Microsoft.Quantum.Canon.ApplyToMostC
    operation ApplyToMostCA<'T>(op : ('T[] => Unit is Adj+Ctl), targets : 'T[]) : Unit is Adj+Ctl {
        op(Most(targets));
    }
}
