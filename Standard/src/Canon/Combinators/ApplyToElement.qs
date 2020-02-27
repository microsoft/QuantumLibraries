// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {

    /// # Summary
    /// Applies an operation to a given element of an array.
    ///
    /// # Description
    /// Given an operation `op`, an index `index`, and an array of targets `targets`,
    /// applies `op(targets[index])`.
    ///
    /// # Input
    /// ## op
    /// An operation to be applied.
    /// ## index
    /// An index into the array of targets.
    /// ## target
    /// An array of possible targets for `op`.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The input type of the operation to be applied.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApplyToElementC
    /// - Microsoft.Quantum.Canon.ApplyToElementA
    /// - Microsoft.Quantum.Canon.ApplyToElementCA
    operation ApplyToElement<'T> (op : ('T => Unit), index : Int, targets : 'T[]) : Unit {
        op(targets[index]);
    }

    /// # Summary
    /// Applies an operation to a given element of an array.
    ///
    /// # Description
    /// Given an operation `op`, an index `index`, and an array of targets `targets`,
    /// applies `op(targets[index])`.
    ///
    /// # Input
    /// ## op
    /// An operation to be applied.
    /// ## index
    /// An index into the array of targets.
    /// ## target
    /// An array of possible targets for `op`.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The input type of the operation to be applied.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApplyToElement
    /// - Microsoft.Quantum.Canon.ApplyToElementA
    /// - Microsoft.Quantum.Canon.ApplyToElementCA
    operation ApplyToElementC<'T> (op : ('T => Unit is Ctl), index : Int, targets : 'T[]) : Unit is Ctl {
        op(targets[index]);
    }

    /// # Summary
    /// Applies an operation to a given element of an array.
    ///
    /// # Description
    /// Given an operation `op`, an index `index`, and an array of targets `targets`,
    /// applies `op(targets[index])`.
    ///
    /// # Input
    /// ## op
    /// An operation to be applied.
    /// ## index
    /// An index into the array of targets.
    /// ## target
    /// An array of possible targets for `op`.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The input type of the operation to be applied.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApplyToElement
    /// - Microsoft.Quantum.Canon.ApplyToElementC
    /// - Microsoft.Quantum.Canon.ApplyToElementCA
    operation ApplyToElementA<'T> (op : ('T => Unit is Adj), index : Int, targets : 'T[]) : Unit is Adj {
        op(targets[index]);
    }

    /// # Summary
    /// Applies an operation to a given element of an array.
    ///
    /// # Description
    /// Given an operation `op`, an index `index`, and an array of targets `targets`,
    /// applies `op(targets[index])`.
    ///
    /// # Input
    /// ## op
    /// An operation to be applied.
    /// ## index
    /// An index into the array of targets.
    /// ## target
    /// An array of possible targets for `op`.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The input type of the operation to be applied.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApplyToElement
    /// - Microsoft.Quantum.Canon.ApplyToElementC
    /// - Microsoft.Quantum.Canon.ApplyToElementA
    operation ApplyToElementCA<'T> (op : ('T => Unit is Adj + Ctl), index : Int, targets : 'T[]) : Unit is Adj + Ctl {
        op(targets[index]);
    }

}
