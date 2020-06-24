namespace Microsoft.Quantum.Math {

    /// # Summary
    /// Represents a rational number of the form `p/q`. Integer `p` is
    /// the first element of the tuple and `q` is the second element
    /// of the tuple.
    newtype Fraction = (
        Numerator: Int,
        Denominator: Int
    );

    /// # Summary
    /// Represents a rational number of the form `p/q`. Integer `p` is
    /// the first element of the tuple and `q` is the second element
    /// of the tuple.
    newtype BigFraction = (
        Numerator: BigInt,
        Denominator: BigInt
    );

}
