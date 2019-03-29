// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.AmplitudeAmplification {
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Primitive;

    /// # Summary
    /// Implementation of <xref:microsoft.quantum.canon.obliviousoraclefromdeterministicstateoracle>.
    operation _ObliviousOracleFromDeterministicStateOracle (ancillaOracle : DeterministicStateOracle, signalOracle : ObliviousOracle, ancillaRegister : Qubit[], systemRegister : Qubit[]) : Unit
    {
        body (...)
        {
            ancillaOracle!(ancillaRegister);
            signalOracle!(ancillaRegister, systemRegister);
        }
        
        adjoint invert;
        controlled distribute;
        controlled adjoint distribute;
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
    /// - Microsoft.Quantum.Canon.DeterministicStateOracle
    /// - Microsoft.Quantum.Canon.ObliviousOracle
    function ObliviousOracleFromDeterministicStateOracle (ancillaOracle : DeterministicStateOracle, signalOracle : ObliviousOracle) : ObliviousOracle
    {
        return ObliviousOracle(_ObliviousOracleFromDeterministicStateOracle(ancillaOracle, signalOracle, _, _));
    }
    
    
    /// # Summary
    /// Implementation of <xref:microsoft.quantum.canon.deterministicstateoraclefromstateoracle>.
    operation _DeterministicStateOracleFromStateOracle (idxFlagQubit : Int, stateOracle : StateOracle, startQubits : Qubit[]) : Unit
    {
        body (...)
        {
            stateOracle!(idxFlagQubit, startQubits);
        }
        
        adjoint invert;
        controlled distribute;
        controlled adjoint distribute;
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
    function DeterministicStateOracleFromStateOracle (idxFlagQubit : Int, stateOracle : StateOracle) : DeterministicStateOracle
    {
        return DeterministicStateOracle(_DeterministicStateOracleFromStateOracle(idxFlagQubit, stateOracle, _));
    }
    
    
    /// # Summary
    /// Implementation of <xref:microsoft.quantum.canon.stateoraclefromdeterministicstateoracle>.
    operation _StateOracleFromDeterministicStateOracle (idxFlagQubit : Int, oracleStateDeterministic : DeterministicStateOracle, qubits : Qubit[]) : Unit
    {
        body (...)
        {
            oracleStateDeterministic!(qubits);
        }
        
        adjoint invert;
        controlled distribute;
        controlled adjoint distribute;
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
    /// Implementation of <xref:microsoft.quantum.canon.reflectionstart>.
    operation _ReflectionStart (phase : Double, qubits : Qubit[]) : Unit
    {
        body (...)
        {
            RAll0(phase, qubits);
        }
        
        adjoint invert;
        controlled distribute;
        controlled adjoint distribute;
    }
    
    
    /// # Summary
    /// Constructs a reflection about the all-zero string |0...0〉, which is the typical input state to amplitude amplification.
    ///
    /// # Output
    /// A `ReflectionOracle` that reflects about the state $\ket{0\cdots 0}$.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ReflectionOracle
    function ReflectionStart () : ReflectionOracle
    {
        return ReflectionOracle(_ReflectionStart(_, _));
    }
    
    
    /// # Summary
    /// Implementation of <xref:microsoft.quantum.canon.reflectionoraclefromdeterministicstateoracle>.
    operation ReflectionOracleFromDeterministicStateOracleImpl (phase : Double, oracle : DeterministicStateOracle, systemRegister : Qubit[]) : Unit
    {
        body (...)
        {
            WithCA(Adjoint oracle!, RAll0(phase, _), systemRegister);
        }
        
        adjoint invert;
        controlled distribute;
        controlled adjoint distribute;
    }
    
    
    /// # Summary
    /// Constructs reflection about a given state from an oracle.
	///
	/// Given the oracle $O$ of type
    /// <xref:microsoft.quantum.canon.deterministicstateoracle>, 
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
    /// - Microsoft.Quantum.Canon.DeterministicStateOracle
    /// - Microsoft.Quantum.Canon.ReflectionOracle
    function ReflectionOracleFromDeterministicStateOracle (oracle : DeterministicStateOracle) : ReflectionOracle
    {
        return ReflectionOracle(ReflectionOracleFromDeterministicStateOracleImpl(_, oracle, _));
    }
    
    
    /// # Summary
    /// Implementation of <xref:microsoft.quantum.canon.targetstatereflectionoracle>.
    operation TargetStateReflectionOracleImpl (phase : Double, idxFlagQubit : Int, qubits : Qubit[]) : Unit {
        body (...) {
            R1(phase, qubits[idxFlagQubit]);
        }

        adjoint invert;
        controlled distribute;
        controlled adjoint distribute;
    }

    /// # Summary
    /// Constructs a `ReflectionOracle` about the target state uniquely marked by the flag qubit.
	///
	/// The target state has a single qubit set to 1, and all others 0: $\ket{1}_f$.
    ///
    /// # Input
    /// ## idxFlagQubit
    /// Index to flag qubit $f$ of oracle.
    ///
    /// # Output
    /// A `ReflectionOracle` that reflects about the state marked by $\ket{1}_f$.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ReflectionOracle
    function TargetStateReflectionOracle (idxFlagQubit : Int) : ReflectionOracle {
        return ReflectionOracle(TargetStateReflectionOracleImpl(_, idxFlagQubit, _));
    }

}


