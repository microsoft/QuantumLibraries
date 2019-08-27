// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Oracles {
    open Microsoft.Quantum.Canon;

    /// # Summary
    /// Implementation of <xref:microsoft.quantum.canon.obliviousoraclefromdeterministicstateoracle>.
    operation _ObliviousOracleFromDeterministicStateOracle (ancillaOracle : DeterministicStateOracle, signalOracle : ObliviousOracle, ancillaRegister : Qubit[], systemRegister : Qubit[]) : Unit
    is Adj + Ctl {
        ancillaOracle::Apply(ancillaRegister);
        signalOracle::Apply(ancillaRegister, systemRegister);
    }


    /// # Summary
    /// Given an oblivious oracle, returns an oblivious oracle that calls a
    /// given determinsitic state oracle on its auxillary register.
    ///
    /// # Input
    /// ## statePreparationOracle
    /// A state preparation oracle.
    /// ## signalOracle
    /// An oblivious oracle.
    ///
    /// # Output
    /// An oblivious oracle that when applied, first applies
    /// `statePreparationOracle` on the auxillary register, and then applies
    /// the original `signalOracle` jointly on the auxillary and system
    /// registers.
    ///
    /// # Remarks
    /// Let `statePreparationOracle` be represented by a unitary operator $A$
    /// acting on the state of an auxillary register $a$, and let `signalOracle`
    /// be represented by a unitary operator $U$ acting jointly on a register
    /// $a, s$. Then the oracle returned by this function can be represented by
    /// a unitary operator $O = U \cdot (A \otimes \id_s)$.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.Bound
    /// - Microsoft.Quantum.Canon.DeterministicStateOracle
    /// - Microsoft.Quantum.Canon.ObliviousOracle
    function BoundWithDeterministicStateOracle(
        statePreparationOracle : DeterministicStateOracle,
        signalOracle : ObliviousOracle
    ) : ObliviousOracle {
        return ObliviousOracle(_ObliviousOracleFromDeterministicStateOracle(statePreparationOracle, signalOracle, _, _));
    }


    /// # Summary
    /// Implementation of <xref:microsoft.quantum.canon.deterministicstateoraclefromstateoracle>.
    operation _DeterministicStateOracleFromStateOracle (idxFlagQubit : Int, stateOracle : StateOracle, startQubits : Qubit[]) : Unit
    is Adj + Ctl {
        stateOracle::Apply(idxFlagQubit, startQubits);
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
    /// - Microsoft.Quantum.Canon.StateOracle
    /// - Microsoft.Quantum.Canon.DeterministicStateOracle
    function DeterministicStateOracleFromStateOracle (idxFlagQubit : Int, stateOracle : StateOracle) : DeterministicStateOracle {
        return DeterministicStateOracle(_DeterministicStateOracleFromStateOracle(idxFlagQubit, stateOracle, _));
    }


    /// # Summary
    /// Implementation of <xref:microsoft.quantum.canon.stateoraclefromdeterministicstateoracle>.
    operation _StateOracleFromDeterministicStateOracle (idxFlagQubit : Int, oracleStateDeterministic : DeterministicStateOracle, qubits : Qubit[]) : Unit
    is Adj + Ctl {
        oracleStateDeterministic::Apply(qubits);
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
    /// - Microsoft.Quantum.Canon.DeterministicStateOracle
    /// - Microsoft.Quantum.Canon.StateOracle
    function StateOracleFromDeterministicStateOracle (deterministicStateOracle : DeterministicStateOracle) : StateOracle
    {
        return StateOracle(_StateOracleFromDeterministicStateOracle(_, deterministicStateOracle, _));
    }


    /// # Summary
    /// Implementation of <xref:microsoft.quantum.canon.reflectionoraclefromdeterministicstateoracle>.
    operation ReflectionOracleFromDeterministicStateOracleImpl (phase : Double, oracle : DeterministicStateOracle, systemRegister : Qubit[]) : Unit
    {
        body (...)
        {
            ApplyWithCA(Adjoint oracle!, RAll0(phase, _), systemRegister);
        }

        adjoint invert;
        controlled distribute;
        controlled adjoint distribute;
    }


    /// # Summary
    /// Constructs reflection about a given state from an oracle.
	///
	/// Given the oracle $O$ of type
    /// <xref:microsoft.quantum.oracles.deterministicstateoracle>,
	/// the result of this function is a reflection around the state $\ket{\psi}$
	/// where $O\ket{0} = \ket{\psi}$.
    ///
    /// # Input
    /// ## oracle
    /// Oracle of type "DeterministicStateOracle"
    ///
    /// # Output
    /// A `ReflectionOracle` that reflects about the state $\ket{\psi}$.
    ///
    /// # See Also
    /// - DeterministicStateOracle
    /// - ReflectionOracle
    function ReflectionOracleFromDeterministicStateOracle (oracle : DeterministicStateOracle) : ReflectionOracle
    {
        return ReflectionOracle(ReflectionOracleFromDeterministicStateOracleImpl(_, oracle, _));
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
        return DiscreteOracle(OperationPowImplCA(blackBoxOracle, _, _));
    }

}
