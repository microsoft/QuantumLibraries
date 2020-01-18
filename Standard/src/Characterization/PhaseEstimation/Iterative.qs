// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Characterization {
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Convert;
    open Microsoft.Quantum.Oracles;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Measurement;

    // Design notes:
    //     The APIs for the iterative and quantum phase estimation algorithms are
    //     parameterized in terms of discrete-time oracles as defined in the Microsoft.Quantum.Oracles
    //     namespace. Constructing such oracle definitions can be done with operations
    //     acting on other operations, making heavy use of partial application internally.

    // E.g.:
    //     let DiscreteOracle = OracleToDiscrete(U);
    //     DiscretePhaseEstimationIteration(oracle, pow, theta, targetState, control);
    //     let datum = M(control);

    //     This design then enables providing more efficient implementations of U^m
    //     to be provided by the user for specific U, while providing a sensible
    //     "default" for operations which cannot be fast-forwarded by taking advantage
    //     of their definitions.

    //     Finally, we also ensure that quantum arguments are placed last to follow
    //     standard conventions for partial application. In particular, this allows for
    //     the discrete and continuous phase estimation iterations to be used in
    //     an allocate-op-measure pattern; we may soon want to abstract that away.

    /// # Summary
    /// Performs a single iteration of an iterative (classically-controlled) phase
    /// estimation algorithm using integer powers of a unitary oracle.
    ///
    /// # Input
    /// ## oracle
    /// Operation acting on an integer and a register,
    /// such that $U^m$ is applied to the given register, where $U$ is the unitary
    /// whose phase is to be estimated, and where $m$ is the integer power
    /// given to the oracle
    /// ## power
    /// Number of times to apply the given unitary oracle.
    /// ## theta
    /// Angle by which to invert the phase on the control qubit before
    /// acting on the target state.
    /// ## targetState
    /// Register of state acted upon by the given unitary oracle.
    /// ## controlQubit
    /// An auxillary qubit used to control the application of the given oracle.
    operation DiscretePhaseEstimationIteration(
        oracle : DiscreteOracle,
        power : Int,
        theta : Double,
        targetState : Qubit[],
        controlQubit : Qubit
    )
    : Unit is Adj + Ctl {
        // NB: We accept the control qubit as input so that we can allow for this operation
        //     to subject to the adjoint and control modifiers (that is, such that we do not need
        //     a return statement, but rather *act* on the given qubits).
        // Find the actual inversion angle by rescaling with the power of the
        // oracle.
        let inversionAngle = -theta * IntAsDouble(power);

        // Prepare the control qubit, using the within/apply block to
        // return the control qubit to the appropriate measurement basis.
        within {
            H(controlQubit);
        } apply {
            Rz(inversionAngle, controlQubit);
            Controlled oracle!([controlQubit], (power, targetState));
        }
    }


    /// # Summary
    /// Performs a single iteration of an iterative (classically-controlled) phase
    /// estimation algorithm using arbitrary real powers of a unitary oracle.
    ///
    /// # Input
    /// ## oracle
    /// Operation acting on a double $t$ and a register,
    /// such that $U^t$ is applied to the given register, where $U$ is the unitary
    /// whose phase is to be estimated, and where $t$ is the integer power
    /// given to the oracle
    /// ## time
    /// Number of times to apply the given unitary oracle.
    /// ## theta
    /// Angle by which to invert the phase on the control qubit before
    /// acting on the target state.
    /// ## targetState
    /// Register of state acted upon by the given unitary oracle.
    operation ContinuousPhaseEstimationIteration(
        oracle : ContinuousOracle,
        time : Double,
        theta : Double,
        targetState : Qubit[],
        controlQubit : Qubit
    ) : Unit is Adj + Ctl
    {
        let inversionAngle = -(theta * time);

        // Prepare the control qubit, using the within/apply block to
        // return the control qubit to the appropriate measurement basis.
        within {
            H(controlQubit);
        } apply {
            Rz(inversionAngle, controlQubit);
            Controlled oracle!([controlQubit], (time, targetState));
        }
    }

}


