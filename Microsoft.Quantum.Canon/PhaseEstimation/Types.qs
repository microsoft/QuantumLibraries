// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {

    /// # Summary
    /// Represents a discrete-time oracle $U^m$ for a fixed operation $U$
    /// and a non-negative integer $m$.
    newtype DiscreteOracle = ((Int, Qubit[]) => ():Adjoint,Controlled);

    /// # Summary
    /// Represents a continuous-time oracle
    /// $U(\delta t) : \ket{\psi(t)} \mapsto \ket{\psi(t + \delta t)}
    /// for all times $t$, where $U$ is a fixed operation, and where
    /// and $\delta t$ is a non-negative real number.
    newtype ContinuousOracle = ((Double, Qubit[]) => ():Adjoint,Controlled);

    // FIXME: Need a better name for this.
    /// # Summary
    /// Given an operation representing a "black-box" oracle, implements
    /// a discrete-time oracle by repeating the given oracle multiple times.
    ///
    /// # Input
    /// ## blackBoxOracle
    /// "Black-box" oracle to perform powers of as a discrete-time oracle.
    ///
    /// # Example
    /// `OracleToDiscrete(U)(3, target)` is equivalent to `U(target)`
    /// repeated three times.
    operation OracleToDiscrete(blackBoxOracle : (Qubit[] => (): Adjoint, Controlled))  : DiscreteOracle
    {
        body {
            let oracle = DiscreteOracle(OperationPowImplAC(blackBoxOracle, _, _));
            return oracle;
        }
    }

}
