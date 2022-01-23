// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Oracles {
    open Microsoft.Quantum.Canon;

    /// # Summary
    /// Implementation of <xref:Microsoft.Quantum.Oracles.ObliviousOracleFromDeterministicStateOracle>.
    internal operation ApplyObliviousOracleFromDeterministicStateOracle(ancillaOracle : DeterministicStateOracle, signalOracle : ObliviousOracle, ancillaRegister : Qubit[], systemRegister : Qubit[])
    : Unit is Adj + Ctl {
        ancillaOracle!(ancillaRegister);
        signalOracle!(ancillaRegister, systemRegister);
    }


    /// # Summary
    /// Combines the oracles `DeterministicStateOracle` and `ObliviousOracle`.
    ///
    /// # Input
    /// ## ancillaOracle
    /// A state preparation oracle $A$ of type `DeterministicStateOracle` acting on register $a$.
    /// ## signalOracle
    /// A oracle $U$ of type `ObliviousOracle` acting jointly on register $a,s$.
    ///
    /// # Output
    /// An oracle $O=UA$ of type `ObliviousOracle`.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Oracles.DeterministicStateOracle
    /// - Microsoft.Quantum.Oracles.ObliviousOracle
    function ObliviousOracleFromDeterministicStateOracle(ancillaOracle : DeterministicStateOracle, signalOracle : ObliviousOracle)
    : ObliviousOracle {
        return ObliviousOracle(ApplyObliviousOracleFromDeterministicStateOracle(ancillaOracle, signalOracle, _, _));
    }


    /// # Summary
    /// Implementation of <xref:Microsoft.Quantum.Oracles.DeterministicStateoracleFromStateOracle>.
    internal operation ApplyDeterministicStateOracleFromStateOracle (idxFlagQubit : Int, stateOracle : StateOracle, startQubits : Qubit[])
    : Unit is Adj + Ctl {
        stateOracle!(idxFlagQubit, startQubits);
    }


    /// # Summary
    /// Converts an oracle of type `StateOracle` to `DeterministicStateOracle`.
    ///
    /// # Input
    /// ## idxFlagQubit
    /// The index to the flag qubit of the `stateOracle` $A$,
    /// which explicitly acts on two registers: the flag $f$ and the system
    /// $s$, e.g. $A\ket{0}\_f\ket{\psi}\_s$.
    /// ## stateOracle
    /// A state preparation oracle $A$ of type `StateOracle`.
    ///
    /// # Output
    /// The same state preparation oracle $A$, but now of type
    /// `DeterministicStateOracle`, so it acts on a register where $a,s$ no
    /// longer explicitly separate, e.g.  $A\ket{0\psi}\_{as}$.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Oracles.StateOracle
    /// - Microsoft.Quantum.Oracles.DeterministicStateOracle
    function DeterministicStateOracleFromStateOracle (idxFlagQubit : Int, stateOracle : StateOracle)
    : DeterministicStateOracle {
        return DeterministicStateOracle(ApplyDeterministicStateOracleFromStateOracle(idxFlagQubit, stateOracle, _));
    }


    /// # Summary
    /// Implementation of <xref:Microsoft.Quantum.Canon.StateOracleFromDeterministicStateOracle>.
    internal operation ApplyStateOracleFromDeterministicStateOracle(idxFlagQubit : Int, oracleStateDeterministic : DeterministicStateOracle, qubits : Qubit[])
    : Unit is Adj + Ctl {
        oracleStateDeterministic!(qubits);
    }


    /// # Summary
    /// Converts an oracle of type `DeterministicStateOracle` to `StateOracle`.
    ///
    /// # Input
    /// ## deterministicStateOracle
    /// A state preparation oracle $A$ of type `DeterministicStateOracle`.
    ///
    /// # Output
    /// The same state preparation oracle $A$, but now of type
    /// `StateOracle`. Note that the flag index in this `StateOracle` is a
    /// dummy variable and has no effect.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Oracles.DeterministicStateOracle
    /// - Microsoft.Quantum.Oracles.StateOracle
    function StateOracleFromDeterministicStateOracle (deterministicStateOracle : DeterministicStateOracle) : StateOracle
    {
        return StateOracle(ApplyStateOracleFromDeterministicStateOracle(_, deterministicStateOracle, _));
    }


    /// # Summary
    /// Implementation of <xref:Microsoft.Quantum.Canon.ReflectionOracleFromDeterministicStateOracle>.
    internal operation _ReflectionOracleFromDeterministicStateOracle(phase : Double, oracle : DeterministicStateOracle, systemRegister : Qubit[])
    : Unit is Adj + Ctl {
        ApplyWithCA(Adjoint oracle!, RAll0(phase, _), systemRegister);
    }

    /// # Summary
    /// Constructs reflection about a given state from an oracle.
    ///
    /// # Description
    /// Given a deterministic state preparation oracle represented by a unitary
    /// matrix $O$,
    /// the result of this function is an oracle that applies a reflection
    /// around the state $\ket{\psi}$ prepared by the oracle $O$; that is,
    /// the state $\ket{\psi}$ such that $O\ket{0} = \ket{\psi}$.
    ///
    /// # Input
    /// ## oracle
    /// An oracle that prepares copies of the state $\ket{\psi}$.
    ///
    /// # Output
    /// An oracle that reflects about the state $\ket{\psi}$.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Oracles.DeterministicStateOracle
    /// - Microsoft.Quantum.Oracles.ReflectionOracle
    function ReflectionOracleFromDeterministicStateOracle(oracle : DeterministicStateOracle)
    : ReflectionOracle {
        return ReflectionOracle(_ReflectionOracleFromDeterministicStateOracle(_, oracle, _));
    }

    /// # Summary
    /// Given an operation representing a "black-box" oracle, returns a discrete-time oracle
    /// which represents the "black-box" oracle repeated multiple times.
    ///
    /// # Input
    /// ## blackBoxOracle
    /// The operation to be exponentiated.
    ///
    /// # Output
    /// An operation partially applied over the "black-box" oracle representing the discrete-time oracle
    ///
    /// # Remarks
    /// ## Example
    /// `OracleToDiscrete(U)(3, target)` is equivalent to `U(target)` repeated three times.
    function OracleToDiscrete (blackBoxOracle : (Qubit[] => Unit is Adj + Ctl)) : DiscreteOracle {
        return DiscreteOracle(ApplyOperationRepeatedlyCA(blackBoxOracle, _, _));
    }

}
