// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Logical {
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Canon;

    /// # Summary
    /// Returns the Boolean negation of a value.
    ///
    /// # Input
    /// ## value
    /// The value to be negated.
    ///
    /// # Output
    /// `true` if and only if `value` is `false`.
    ///
    /// # Remarks
    /// The following are equivalent:
    /// ```Q#
    /// let x = not value;
    /// let x = Not(value);
    /// ```
    function Not(value : Bool) : Bool {
        return not value;
    }

    /// # Summary
    /// Returns the Boolean conjunction of two values.
    ///
    /// # Input
    /// ## a
    /// The first value to be considered.
    /// ## b
    /// The second value to be considered.
    ///
    /// # Output
    /// `true` if and only if `a` and `b` are both `true`.
    ///
    /// # Remarks
    /// Unlike the `and` operator, this function does not short-circuit, such that
    /// both inputs are fully evaluated.
    ///
    /// Up to short-circuiting behavior, the following are equivalent:
    /// ```Q#
    /// let x = a and b;
    /// let x = And(a, b);
    /// ```
    function And(a : Bool, b : Bool) : Bool {
        return a and b;
    }

    /// # Summary
    /// Returns the Boolean disjunction of two values.
    ///
    /// # Input
    /// ## a
    /// The first value to be considered.
    /// ## b
    /// The second value to be considered.
    ///
    /// # Output
    /// `true` if and only if either `a` or `b` are `true`.
    ///
    /// # Remarks
    /// Unlike the `or` operator, this function does not short-circuit, such that
    /// both inputs are fully evaluated.
    ///
    /// Up to short-circuiting behavior, the following are equivalent:
    /// ```Q#
    /// let x = a or b;
    /// let x = Or(a, b);
    /// ```
    function Or(a : Bool, b : Bool) : Bool {
        return a or b;
    }

    /// # Summary
    /// Returns the Boolean exclusive disjunction of two values.
    ///
    /// # Input
    /// ## a
    /// The first value to be considered.
    /// ## b
    /// The second value to be considered.
    ///
    /// # Output
    /// `true` if and only if exactly one of `a` and `b` is `true`.
    function Xor(a : Bool, b : Bool) : Bool {
        return (a or b) and ((not a) or (not b));
    }

    /// # Summary
    /// Returns one of two values, depending on the value of a Boolean condition.
    ///
    /// # Input
    /// ## condition
    /// A condition used to control which input is returned.
    /// ## ifTrue
    /// The value to be returned when `condition` is `true`.
    /// ## ifFalse
    /// The value to be returned when `condition` is `false`.
    ///
    /// # Output
    /// `ifTrue` if `condition` is `true`, and `ifFalse` otherwise.
    ///
    /// # Remarks
    /// Unlike the `?|` operator, this function does not short-circuit, such that
    /// both inputs are fully evaluated.
    ///
    /// Up to short-circuiting behavior, the following are equivalent:
    /// ```Q#
    /// let x = condition ? ifTrue | ifFalse;
    /// let x = Conditioned(condition, ifTrue, ifFalse);
    /// ```
    function Conditioned<'T>(condition : Bool, ifTrue : 'T, ifFalse : 'T) : 'T {
        return condition ? ifTrue | ifFalse;
    }

}
