namespace Microsoft.Quantum.Canon {

    operation OperationPowImpl<'T>(oracle : ('T => ()), power : Int, target : 'T)  : ()
    {
        body {
            for (idxApplication in 0..power - 1) {
                oracle(target);
            }
        }
    }

    operation OperationPowImplAC<'T>(oracle : ('T => () : Controlled,Adjoint), power : Int, target : 'T)  : ()
    {
        body {
            for (idxApplication in 0..power - 1) {
                oracle(target);
            }
        }

        adjoint auto
        controlled auto
        controlled adjoint auto
    }

    /// # Summary
    /// Given an operation representing a gate $U$, returns a new operation
    /// $U^m$ for a power $m$.
    ///
    /// ## oracle
    /// An operation $U$ representing the gate to be repeated.
    /// ## power
    /// The number of times that $U$ is to be repeated.
    ///
    /// # Output
    /// A new operation representing $U^m$, where $m = \texttt{power}$.
    function OperationPow<'T>(oracle : ('T => ()), power : Int)  : ('T => ())
    {
        return OperationPowImpl(oracle, power, _);
    }

    // TODO: A, C, AC variants.

}