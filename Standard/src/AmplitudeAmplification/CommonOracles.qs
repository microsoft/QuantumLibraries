﻿// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.AmplitudeAmplification {
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Oracles;

    /// # Summary
    /// Constructs a reflection about the all-zero string |0...0〉, which is the typical input state to amplitude amplification.
    ///
    /// # Output
    /// A `ReflectionOracle` that reflects about the state $\ket{0\cdots 0}$.
    ///
    /// # See Also
    /// - Microsoft.Quantum.AmplitudeAmplification.ReflectionOracle
    function ReflectionStart() : ReflectionOracle {
        return ReflectionOracle(RAll0);
    }

    /// # Summary
    /// Implementation of <xref:Microsoft.Quantum.AmplitudeAmplification.TargetStateReflectionOracle>.
    internal operation ApplyTargetStateReflectionOracle(phase : Double, idxFlagQubit : Int, qubits : Qubit[])
    : Unit is Adj + Ctl {
        R1(phase, qubits[idxFlagQubit]);
    }

    /// # Summary
    /// Constructs a `ReflectionOracle` about the target state uniquely marked by the flag qubit.
    ///
    /// # Description
    /// The oracle returned by this function reflects about the $\ket{1}$ state
    /// of a flag qubit (given by the index `idxFlagQubit`).
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
    function TargetStateReflectionOracle(idxFlagQubit : Int) : ReflectionOracle {
        return ReflectionOracle(ApplyTargetStateReflectionOracle(_, idxFlagQubit, _));
    }

}

