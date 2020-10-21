namespace Microsoft.Quantum.Math {

    /// # Summary
    /// Represents a rational number of the form `p/q`. Integer `p` is
    /// the first element of the tuple and `q` is the second element
    /// of the tuple.
    ///
    /// # Named Items
    /// ## Numerator
    /// Numerator of the fraction.
    /// ## Denominator
    /// Denominator of the fraction/
    newtype Fraction = (
        Numerator: Int,
        Denominator: Int
    );

    /// # Summary
    /// Represents a rational number of the form `p/q`. Integer `p` is
    /// the first element of the tuple and `q` is the second element
    /// of the tuple.
    ///
    /// # Named Items
    /// ## Numerator
    /// Numerator of the fraction.
    /// ## Denominator
    /// Denominator of the fraction/
    newtype BigFraction = (
        Numerator: BigInt,
        Denominator: BigInt
    );

}
