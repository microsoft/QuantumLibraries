// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {
	open Microsoft.Quantum.Primitive;

    /// # Summary
    /// This performs a phase shift operation about the state $\ket{1\cdots 1}\bra{1\cdots 1}$.
    operation RAll1( phase: Double, qubits: Qubit[] ) : ()
    {
	    body {
		    let nQubits = Length(qubits);
		    let flagQubit = qubits[0];
		    let systemRegister = qubits[1..nQubits-1];

		    (Controlled R1(phase, _))(systemRegister, flagQubit);

	    }

	    adjoint auto
	    controlled auto
	    controlled adjoint auto
    }

    /// # Summary
    /// This performs a phase shift operation about the state $\ket{0\cdots 0}\bra{0\cdots 0}$,
    /// $$
    /// \begin{align}
    ///     \ket{0\cdots 0} \mapsto e^{i \theta} \ket{0\cdots 0}.
    /// \end{align}
    /// $$
    ///
    /// # Input
    /// ## phase
    /// The phase $\theta$ assigned to the state $\ket{0\cdots 0}$.
    /// ## qubits
    /// The register whose state is to be rotated.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.RAll1
    operation RAll0( phase: Double, qubits: Qubit[] ) : ()
    {
	    body {

		    WithCA(ApplyToEachCA(X, _), RAll1(phase, _), qubits);

	    }

	    adjoint auto
	    controlled auto
	    controlled adjoint auto
    }


    /// # Summary
    /// Implementation of @"quantum.microsoft.canon.ObliviousOracleFromDeterministicStateOracle".
    operation _ObliviousOracleFromDeterministicStateOracle(ancillaOracle : DeterministicStateOracle, signalOracle : ObliviousOracle, ancillaRegister: Qubit[], systemRegister: Qubit[]) : (){
	    body{
			    ancillaOracle(ancillaRegister);
			    signalOracle(ancillaRegister, systemRegister);
	    }
	    adjoint auto
	    controlled auto
	    adjoint controlled auto
    }
    /// # Summary
    /// Combines the oracles DeterministicStateOracle and ObliviousOracle.
    function ObliviousOracleFromDeterministicStateOracle(ancillaOracle : DeterministicStateOracle, signalOracle : ObliviousOracle) : ObliviousOracle{
	    return ObliviousOracle(_ObliviousOracleFromDeterministicStateOracle(ancillaOracle, signalOracle,_,_));
    }

    /// # Summary
    /// Implementation of @"quantum.microsoft.canon.DeterministicStateOracleFromStateOracle".
    operation _DeterministicStateOracleFromStateOracle(idxFlagQubit: Int, stateOracle : StateOracle, startQubits: Qubit[]) : (){
	    body{
		    stateOracle(idxFlagQubit, startQubits);
	    }
	    adjoint auto
	    controlled auto
	    adjoint controlled auto
    }
    
    /// # Summary
    /// Converts an oracle of type StateOracle to DeterministicStateOracle.
    function DeterministicStateOracleFromStateOracle(idxFlagQubit: Int, stateOracle : StateOracle) : DeterministicStateOracle{
	    return DeterministicStateOracle(_DeterministicStateOracleFromStateOracle(idxFlagQubit, stateOracle,_));
    }

    /// # Summary
    /// Implementation of @"quantum.microsoft.canon.StateOracleFromDeterministicStateOracle".
    operation _StateOracleFromDeterministicStateOracle(idxFlagQubit : Int, oracleStateDeterministic : DeterministicStateOracle, qubits: Qubit[]): ()
    {
	    body {
		    oracleStateDeterministic(qubits);
	    }
	    adjoint auto
	    controlled auto
	    controlled adjoint auto
    }

    /// # Summary
    /// Converts an oracle of type DeterministicStateOracle to StateOracle.    
    function StateOracleFromDeterministicStateOracle(oracleStateDeterministic : DeterministicStateOracle) : StateOracle {
	    return StateOracle(_StateOracleFromDeterministicStateOracle(_, oracleStateDeterministic,_));
    }

    /// # Summary
    /// Implementation of @"quantum.microsoft.canon.ReflectionStart".
    operation _ReflectionStart(phase: Double, qubits: Qubit[] ) : () {
	    body {
		    RAll0(phase, qubits );
	    }
	    adjoint auto
	    controlled auto
	    adjoint controlled auto
    }

    /// # Summary
    /// Constructs a reflection about the all-zero string $\ket{0\cdots 0}$, which is the typical input state to amplitude amplification.
    ///
    /// # Output
    /// A `ReflectionOracle` that reflects about the state $\ket{0\cdots 0}$.
    function ReflectionStart() : ReflectionOracle {
	    return ReflectionOracle(_ReflectionStart( _, _ ));
    }

    /// # Summary
    /// Implementation of @"microsoft.quantum.canon.ReflectionOracleFromDeterministicStateOracle".
    operation ReflectionOracleFromDeterministicStateOracleImpl(phase: Double, oracle: DeterministicStateOracle, systemRegister: Qubit[]): ()
    {
	    body {
		    WithCA((Adjoint oracle), RAll0(phase, _), systemRegister);
	    }
	    adjoint auto
	    controlled auto
	    adjoint controlled auto
    }

    /// # Summary
    /// Constructs reflection about a some state $\ket{\psi}$ from the oracle of type
    /// @"microsoft.quantum.canon.deterministicstateoracle", where $O\ket{0} = \ket{\psi}$.
    ///
    /// # Input
    /// ## oracle
    /// Oracle of type "DeterministicStateOracle"
    ///
    /// # Output
    /// A `ReflectionOracle` that reflects about the state $\ket{\psi}$.
    function ReflectionOracleFromDeterministicStateOracle(oracle: DeterministicStateOracle): ReflectionOracle
    {
	    return ReflectionOracle(ReflectionOracleFromDeterministicStateOracleImpl(_, oracle, _ ));
    }

    /// # Summary
    /// Implementation of @"microsoft.quantum.canon.TargetStateReflectionOracle".
    operation TargetStateReflectionOracleImpl(phase: Double, idxFlagQubit : Int, qubits: Qubit[]): ()
    {
	    body {
		    R1(phase, qubits[idxFlagQubit]);
	    }

	    adjoint auto
	    controlled auto
	    adjoint controlled auto
    }

    /// # Summary
    /// Constructs reflection about the target state uniquely marked by the flag qubit state $\ket{1}_f$, prepared the oracle of type "StateOracle".
    ///
    /// # Input
    /// ## idxFlagQubit
    /// Index to flag qubit $f$ of oracle.
    ///
    /// # Output
    /// A `ReflectionOracle` that reflects about the state marked by $\ket{1}_f$.
    function TargetStateReflectionOracle(idxFlagQubit : Int): ReflectionOracle
    {
	    return ReflectionOracle(TargetStateReflectionOracleImpl( _ , idxFlagQubit , _ ));
    }

}
