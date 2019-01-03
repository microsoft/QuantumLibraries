// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon
{
    
    /// # Summary
    /// Represents a discrete-time oracle.
	///
	/// This is an oracle that implements $U^m$ for a fixed operation $U$
    /// and a non-negative integer $m$.
    newtype DiscreteOracle = ((Int, Qubit[]) => Unit : Adjoint, Controlled);
    
    /// # Summary
    /// Represents a continuous-time oracle.
	///
	/// This is an oracle that implements 
    /// $U(\delta t) : \ket{\psi(t)} \mapsto \ket{\psi(t + \delta t)}
    /// for all times $t$, where $U$ is a fixed operation, and where
    /// $\delta t$ is a non-negative real number.
    newtype ContinuousOracle = ((Double, Qubit[]) => Unit : Adjoint, Controlled);
    
    
    /// # Summary
    /// Given an operation representing a "black-box" oracle, returns a discrete-time oracle
    /// which represents the "black-box" oracle repeated multiple times.
    ///
    /// # Input
    /// ## blackBoxOracle
    /// The operation to be exponentiated
    ///
    /// # Output
    /// An operation partially applied over the "black-box" oracle representing the discrete-time oracle
    ///
    /// # Remarks
    /// ## Example
    /// `OracleToDiscrete(U)(3, target)` is equivalent to `U(target)` repeated three times.
    operation OracleToDiscrete (blackBoxOracle : (Qubit[] => Unit : Adjoint, Controlled)) : DiscreteOracle
    {
        let oracle = DiscreteOracle(OperationPowImplCA(blackBoxOracle, _, _));
        return oracle;
    }
    
}


