// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.AmplitudeAmplification {
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Convert;
    open Microsoft.Quantum.Math;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Oracles;

    /// # Summary
    /// Oblivious amplitude amplification by specifying partial reflections.
    ///
    /// # Input
    /// ## phases
    /// Phases of partial reflections
    /// ## startStateReflection
    /// Reflection operator about start state of auxiliary register
    /// ## targetStateReflection
    /// Reflection operator about target state of auxiliary register
    /// ## signalOracle
    /// Unitary oracle $O$ of type `ObliviousOracle` that acts jointly on the
    /// auxiliary and system registers.
    /// ## auxiliaryRegister
    /// Auxiliary register
    /// ## systemRegister
    /// System register
    ///
    /// # Remarks
    /// Given a particular auxiliary start state $\ket{\text{start}}\_a$, a
    /// particular auxiliary target state $\ket{\text{target}}\_a$, and any
    /// system state $\ket{\psi}\_s$, suppose that
    /// \begin{align}
    /// O\ket{\text{start}}\_a\ket{\psi}\_s= \lambda\ket{\text{target}}\_a U \ket{\psi}\_s + \sqrt{1-|\lambda|^2}\ket{\text{target}^\perp}\_a\cdots
    /// \end{align}
    /// for some unitary $U$.
    /// By a sequence of reflections about the start and target states on the
    /// auxiliary register interleaved by applications of `signalOracle` and its
    /// adjoint, the success probability of applying U may be altered.
    ///
    /// In most cases, `auxiliaryRegister` is initialized in the state $\ket{\text{start}}\_a$.
    ///
    /// # References
    /// See
    /// - [ *D.W. Berry, A.M. Childs, R. Cleve, R. Kothari, R.D. Somma* ](https://arxiv.org/abs/1312.1414)
    /// for the standard version.
    /// See
    /// - [ *G.H. Low, I.L. Chuang* ](https://arxiv.org/abs/1610.06546)
    /// for a generalization to partial reflections.
    operation ApplyObliviousAmplitudeAmplification(
        phases : ReflectionPhases,
        startStateReflection : ReflectionOracle,
        targetStateReflection : ReflectionOracle,
        signalOracle : ObliviousOracle,
        auxiliaryRegister : Qubit[],
        systemRegister : Qubit[]
    )
    : Unit is Adj + Ctl {
        for ((startPhase, targetPhase) in Zipped(phases!)) {
            if (startPhase != 0.0) {
                startStateReflection::ApplyReflection(
                    startPhase, auxiliaryRegister
                );
            }

            within {
                signalOracle!(auxiliaryRegister, systemRegister);
            } apply {
                if (targetPhase != 0.0) {
                    targetStateReflection!(targetPhase, auxiliaryRegister);
                }
            }
        }
        // This gives us one extra application of Adjoint signalOracle!, so we
        // apply the forward direction at the end.
        signalOracle!(auxiliaryRegister, systemRegister);
    }


    // NB [STYLE]: The name of this operation uses "From" as it is not a type
    //             conversion function ("As"), but something that constructs
    //             an operation from given information in a deterministic
    //             fashion.
    /// # Summary
    /// Returns a unitary that implements oblivious amplitude amplification by specifying for partial reflections.
    function ObliviousAmplitudeAmplificationFromPartialReflections(
        phases : ReflectionPhases,
        startStateReflection : ReflectionOracle,
        targetStateReflection : ReflectionOracle,
        signalOracle : ObliviousOracle
    )
    : ((Qubit[], Qubit[]) => Unit is Adj + Ctl) {
        return ApplyObliviousAmplitudeAmplification(
            phases, startStateReflection, targetStateReflection, signalOracle,
            _, _
        );
    }


    /// # Summary
    /// Oblivious amplitude amplification by oracles for partial reflections.
    ///
    /// # Input
    /// ## phases
    /// Phases of partial reflections
    /// ## startStateOracle
    /// Unitary oracle $A$ that prepares auxiliary start state
    /// ## signalOracle
    /// Unitary oracle $O$ of type `ObliviousOracle` that acts jointly on the
    /// auxiliary and system register
    /// ## idxFlagQubit
    /// Index to single-qubit flag register
    ///
    /// # Output
    /// An operation that implements oblivious amplitude amplification based on partial reflections.
    ///
    /// # Remarks
    /// This imposes stricter conditions on form of the auxiliary start and target states than in `AmpAmpObliviousByReflectionPhases`.
    /// It is assumed that $A\ket{0}\_f\ket{0}\_a= \ket{\text{start}}\_{fa}$ prepares the auxiliary start state $\ket{\text{start}}\_{fa}$ from the computational basis $\ket{0}\_f\ket{0}$.
    /// It is assumed that the target state is marked by $\ket{1}\_f$.
    /// It is assumed that
    /// \begin{align}
    /// O\ket{\text{start}}\_{fa}\ket{\psi}\_s= \lambda\ket{1}\_f\ket{\text{anything}}\_a\ket{\text{target}}\_s U \ket{\psi}\_s + \sqrt{1-|\lambda|^2}\ket{0}\_f\cdots,
    /// \end{align}
    /// for some unitary $U$.
    function ObliviousAmplitudeAmplificationFromStatePreparation(
        phases : ReflectionPhases,
        startStateOracle : DeterministicStateOracle,
        signalOracle : ObliviousOracle,
        idxFlagQubit : Int
    )
    : ((Qubit[], Qubit[]) => Unit is Adj + Ctl) {
        let startStateReflection = ReflectionStart();
        let targetStateReflection = TargetStateReflectionOracle(idxFlagQubit);
        let obliviousSignalOracle = ObliviousOracleFromDeterministicStateOracle(
            startStateOracle, signalOracle
        );
        return ObliviousAmplitudeAmplificationFromPartialReflections(
            phases, startStateReflection, targetStateReflection, obliviousSignalOracle
        );
    }

    /// # Summary
    /// Applies amplitude amplification on a given register, using a given set
    /// of phases and oracles to reflect about the initial and final states.
    ///
    /// # Input
    /// ## phases
    /// A set of phases describing the partial reflections at each step of the
    /// amplitude amplification algorithm. See
    /// @"microsoft.quantum.amplitudeamplification.standardreflectionphases"
    /// for an example.
    /// ## startStateReflection
    /// An oracle that reflects about the initial state.
    /// ## targetStateReflection
    /// An oracle that reflects about the desired final state.
    /// ## target
    /// A register to perform amplitude amplification on.
    operation ApplyAmplitudeAmplification(
        phases : ReflectionPhases,
        startStateReflection : ReflectionOracle,
        targetStateReflection : ReflectionOracle,
        target : Qubit[]
    )
    : Unit is Adj + Ctl {
        // Pass empty qubit array using fact that NoOp does nothing.
        let systemRegister = new Qubit[0];
        let signalOracle = ObliviousOracle(NoOp<(Qubit[], Qubit[])>);
        let op = ObliviousAmplitudeAmplificationFromPartialReflections(
            phases, startStateReflection, targetStateReflection, signalOracle
        );

        op(target, systemRegister);
    }


    /// # Summary
    /// Amplitude amplification by partial reflections.
    ///
    /// # Input
    /// ## phases
    /// Phases of partial reflections
    /// ## startStateReflection
    /// Reflection operator about start state
    /// ## targetStateReflection
    /// Reflection operator about target state
    ///
    /// # Output
    /// An operation that implements amplitude amplification by partial reflections.
    ///
    /// # Remarks
    /// Amplitude amplification is a special case of oblivious amplitude amplification where there are no system qubits and the oblivious oracle is set to identity.
    /// In most cases, `startQubits` is initialized in the state $\ket{\text{start}}\_1$, which is the $-1$ eigenstate of `startStateReflection`.
    function AmplitudeAmplificationFromPartialReflections(
        phases : ReflectionPhases,
        startStateReflection : ReflectionOracle,
        targetStateReflection : ReflectionOracle
    )
    : (Qubit[] => Unit is Adj + Ctl) {
        // Pass empty qubit array using fact that NoOp does nothing.
        let qubitEmpty = new Qubit[0];
        let signalOracle = ObliviousOracle(NoOp<(Qubit[], Qubit[])>);
        return (ObliviousAmplitudeAmplificationFromPartialReflections(
            phases, startStateReflection, targetStateReflection, signalOracle
        ))(_, qubitEmpty);
    }


    // NB [STYLE]: The name of this operation uses "From" as it is not a type
    //             conversion function ("As"), but something that constructs
    //             an operation from given information in a deterministic
    //             fashion.
    /// # Summary
    /// Amplitude amplification by oracles for partial reflections.
    ///
    /// # Input
    /// ## phases
    /// Phases of partial reflections
    /// ## stateOracle
    /// Unitary oracle $A$ that prepares start state
    /// ## idxFlagQubit
    /// Index to flag qubit
    ///
    /// # Output
    /// An operation that implements amplitude amplification by oracles that are
    /// implemented by partial reflections.
    ///
    /// # Remarks
    /// This imposes stricter conditions on form of the start and target states than in `AmpAmpByReflectionPhases`.
    /// It is assumed that the target state is marked by $\ket{1}\_f$.
    /// It is assumed that
    /// \begin{align}
    /// A\ket{0}\_{f}\ket{0}\_s= \lambda\ket{1}\_f\ket{\text{target}}\_s + \sqrt{1-|\lambda|^2}\ket{0}\_f\cdots,
    /// \end{align}
    /// In most cases, `flagQubit` and `auxiliaryRegister` are initialized in the state $\ket{0}\_{f}\ket{0}\_s$.
    function AmplitudeAmplificationFromStatePreparation(
        phases : ReflectionPhases,
        stateOracle : StateOracle,
        idxFlagQubit : Int
    )
    : (Qubit[] => Unit is Adj + Ctl) {
        let systemRegister = new Qubit[0];
        let signalOracle = ObliviousOracle(NoOp<(Qubit[], Qubit[])>);
        let startStateOracle = DeterministicStateOracleFromStateOracle(idxFlagQubit, stateOracle);
        return (ObliviousAmplitudeAmplificationFromStatePreparation(
            phases, startStateOracle, signalOracle, idxFlagQubit
        ))(_, systemRegister);
    }


    /// # Summary
    /// Standard Amplitude Amplification algorithm
    ///
    /// # Input
    /// ## nIterations
    /// Number of iterations $n$ of amplitude amplification
    /// ## stateOracle
    /// Unitary oracle $A$ that prepares start state
    /// ## idxFlagQubit
    /// Index to flag qubit
    ///
    /// # Output
    /// An operation that implements the standard amplitude amplification quantum algorithm
    ///
    /// # Remarks
    /// This is the standard amplitude amplification algorithm obtained by a choice of reflection phases computed by `AmpAmpPhasesStandard`
    /// Assuming that
    /// \begin{align}
    /// A\ket{0}\_{f}\ket{0}\_s= \lambda\ket{1}\_f\ket{\text{target}}\_s + \sqrt{1-|\lambda|^2}\ket{0}\_f\cdots,
    /// \end{align}
    /// this operation prepares the state
    /// \begin{align}
    /// \operatorname{AmpAmpByOracle}\ket{0}\_{f}\ket{0}\_s= \sin((2n+1)\sin^{-1}(\lambda))\ket{1}\_f\ket{\text{target}}\_s + \cdots\ket{0}\_f
    /// \end{align}
    /// In most cases, `flagQubit` and `auxiliaryRegister` is initialized in the state $\ket{0}\_f\ket{0}\_a$.
    ///
    /// # References
    /// - [ *G. Brassard, P. Hoyer, M. Mosca, A. Tapp* ](https://arxiv.org/abs/quant-ph/0005055)
    function StandardAmplitudeAmplification(
        nIterations : Int,
        stateOracle : StateOracle,
        idxFlagQubit : Int
    )
    : (Qubit[] => Unit is Adj + Ctl) {
        let phases = StandardReflectionPhases(nIterations);
        return AmplitudeAmplificationFromStatePreparation(phases, stateOracle, idxFlagQubit);
    }


    /// # Summary
    /// Fixed-Point Amplitude Amplification algorithm
    ///
    /// # Input
    /// ## statePrepOracle
    /// Unitary oracle that prepares the start state.
    /// ## startQubits
    /// Qubit register
    ///
    /// # Remarks
    /// The startQubits must be in the $\ket{0 \cdots 0}$ state. This operation iterates over a number of queries in powers of $2$ until either a maximal number of queries
    /// is reached, or the target state is found.
    operation ApplyFixedPointAmplification(statePrepOracle : StateOracle, startQubits : Qubit[])
    : Unit {
        // Should be a power of 2
        let queriesMax = 999;
        let successMin = 0.99;
        mutable finished = Zero;
        mutable exponentMax = 0;
        mutable exponentCurrent = 0;

        //Complexity: Let \theta = \mathcal{O}(\sqrt{lambda})
        // Number of Measurements = O( Log^2(1/\theta) )
        // Number of Queries = O(1/\theta)
        using (flagQubit = Qubit[1]) {
            let qubits = flagQubit + startQubits;
            let idxFlagQubit = 0;

            repeat {
                if (2 ^ exponentMax > queriesMax) {
                    fail $"Target state not found. Maximum number of queries exceeded.";
                }

                repeat {
                    let queries = 2 ^ exponentCurrent;
                    let phases = FixedPointReflectionPhases(queries, successMin);
                    (AmplitudeAmplificationFromStatePreparation(phases, statePrepOracle, idxFlagQubit))(qubits);
                    set finished = M(flagQubit[0]);
                    set exponentCurrent = exponentCurrent + 1;
                }
                until (finished == One or exponentCurrent > exponentMax)
                fixup {
                    // flagQubit is already in Zero for fixup to apply
                    ResetAll(startQubits);
                }

                set exponentCurrent = 0;
                set exponentMax = exponentMax + 1;
            }
            until (finished == One)
            fixup {
                ResetAll(startQubits);
            }
        }
    }

}


