namespace Microsoft.Quantum.Math {

    /// # Summary
    /// Represents a rational number of the form `p/q`. Integer `p` is
    /// the first element of the tuple and `q` is the second element
    /// of the tuple.
    newtype Fraction = (Int, Int);

    /// # Summary
    /// Represents a rational number of the form `p/q`. Integer `p` is
    /// the first element of the tuple and `q` is the second element
    /// of the tuple.
    newtype BigFraction = (BigInt, BigInt);

}
