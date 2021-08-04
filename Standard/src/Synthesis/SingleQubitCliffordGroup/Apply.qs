// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Synthesis {

    open Microsoft.Quantum.Convert;
    open Microsoft.Quantum.Math;
    open Microsoft.Quantum.Intrinsic;

    // Design notes:
    // Thinking in terms of bounded polymorphism, this would be part of a
    // concept like ApplicableTo<'TInput>,
    //
    // concept 'TOperation is ApplicableToCA<'TInput> when {
    //     operation Apply(op : 'TOperation, input : 'TInput) : Unit is Adj + Ctl;
    // }
    // example <'TInput> ('TInput => Unit is Adj + Ctl) is ApplicableToCA<'TInput> {
    //     operation Apply(op : ('TInput => Unit is Adj + Ctl), input : 'TInput) : Unit is Adj + Ctl {
    //         op(input)
    //     }
    // }
    //
    // operation ApplyToEachCA<'TOperation, 'TInput where 'TOperation is ApplicableToCA<'TInput>(
    //     operation : 'TOperation,
    //     inputs : 'TInput[]
    // )
    // : Unit is Adj + Ctl {
    //     for input in inputs { Apply(op, input); }
    // }
    //
    // Taking that concept as inspiration, we can design the signature here
    // to match, up to the use of a type suffix to denote the particular
    // instance of the ApplicableToCA concept.
    //
    // We will use the "1C" suffix for now.

    internal operation ApplyDirectly1C(op : SingleQubitClifford, target : Qubit) : Unit is Adj + Ctl {
        let cOp = CanonicalForm1C(op);
        // Following this presentation of the Clifford group, we need to apply
        // ω^𝑙 𝐸^𝑖 𝑋^𝑗 𝑆^𝑘, where (𝑖, 𝑗, 𝑘, 𝑙) is the given member of the Clifford
        // group.

        // First, either apply 𝑋 or not depending on 𝑗.
        if   (cOp::X == 1) { X(target); }

        // Apply 𝑆^𝑘 next. We use that 𝑆² = 𝑍 to simplify what operations we
        // need to apply.
        if   (cOp::S == 1) { S(target); }
        elif (cOp::S == 2) { Z(target); }
        elif (cOp::S == 3) { Adjoint S(target); }

        // Next, apply 𝐸^𝑖. Recall that 𝐸 ≔ 𝐻𝑆⁺ω³. Using that, we'll ignore
        // the ω part for now and will fold it in below when we correct the
        // global phase.
        // NOTE: We could optimize this further by collapsing the Adjoint S
        //       here with the S from above in some cases, but we'll leave
        //       that for the compiler.
        for _ in 1..cOp::E { Adjoint S(target); H(target); }

        // Finally, correct the global phase. While this doesn't matter for
        // applying the operation in an uncontrolled fashion, that phase
        // becomes local under the Controlled functor.
        // To do so, we also need to include the global phase we ignored above
        // when we applied the 𝐸 part of this operation.
        let actualOmega = (cOp::Omega + 3 * cOp::E) % 8;

        R(PauliI, 2.0 * PI() * IntAsDouble(actualOmega) / 8.0, target);

    }

    /// # Summary
    /// Given a single-qubit Clifford operator, applies the corresponding operation
    /// to a single qubit.
    ///
    /// # Input
    /// ## op
    /// The Clifford operator to be applied.
    /// ## target
    /// The qubit to which `op` is to be applied as an operation.
    ///
    /// # Example
    /// The following are equivalent, as $ES\omega^5 = (HS^3\omega^3)S\omega^5
    /// = HS^4\omega^8 = H$:
    /// ```qsharp
    /// Apply1C(Identity1C() w/ E <- 1 w/ S <- 1 w/ Omega <- 5), q);
    /// ```
    /// and
    /// ```qsharp
    /// H(q);
    /// ```
    operation Apply1C(op : SingleQubitClifford, target : Qubit) : Unit is Adj + Ctl {
        body (...) {
            ApplyDirectly1C(op, target);
        }
        adjoint (...) {
            ApplyDirectly1C(Inverse1C(op), target);
        }
    }

}
