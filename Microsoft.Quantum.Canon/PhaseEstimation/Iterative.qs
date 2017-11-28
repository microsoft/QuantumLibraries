// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {
    open Microsoft.Quantum.Primitive;

    // Design notes:
    //     The APIs for the iterative and quantum phase estimation algorithms are
    //     parameterized in terms of discrete-time oracles as defined in the OracleTypes.qs
    //     source file. Constructing such oracle definitions can be done with operations
    //     acting on other operations, making heavy use of partial application internally.
    //
    // E.g.:
    //     let DiscreteOracle = OracleToDiscrete(U)
    //     DiscretePhaseEstimationIteration(oracle, pow, theta, eigenstate, control)
    //     let Result datum = Measure(control)
    //
    //     This design then enables providing more efficient implementations of U^m
    //     to be provided by the user for specific U, while providing a sensible
    //     "default" for operations which cannot be fast-forwarded by taking advantage
    //     of their definitions.
    //
    //     Finally, we also ensure that quantum arguments are placed last to follow
    //     standard conventions for partial application. In particular, this allows for
    //     the discrete and continuous phase estimation iterations to be used in
    //     an allocate-op-measure pattern; we may soon want to abstract that away.

    /// Performs a single iteration of an iterative (classically-controlled) phase
    /// estimation algorithm.
    /// <param name="oracle">Operation acting on an integer and a register,
    ///     such that U^m is applied to the given register, where U is the unitary
    ///     whose phase is to be estimated, and where m is the integer power
    ///     given to the oracle.</param>
    /// <param name="eigenstate">Register containing an eigenstate of the given oracle.</param>
    /// <param name="power">Number of times to apply the given unitary oracle.</param>
    /// <param name="theta">Angle by which to invert the phase on the control qubit before
    ///     acting on the eigenstate.</param>
    operation DiscretePhaseEstimationIteration( oracle : DiscreteOracle, power : Int, theta : Double, eigenstate : Qubit[], controlQubit : Qubit)  : ()
    {
        // NB: We accept the control qubit as input so that we can allow for this operation
        //     to subject to the adjoint and control modifiers (that is, such that we do not need
        //     a return statement, but rather *act* on the given qubits).
        body {
            // if (power < 0) {
                // fail 'Oracle power cannot be negative.'
            // }

            // Find the actual inversion angle by rescaling with the power of the
            // oracle.
            let inversionAngle = -theta * Float(power);
            
            // Prepare the control qubit.
            H(controlQubit);
            Rz(inversionAngle, controlQubit);

            // TODO: should this be
            //     controlled(oracle)([controlQubit], power, eigenstate),
            // or
            //     controlled(oracle(power, _))([controlQubit], eigenstate)?
            (Controlled oracle)([controlQubit], (power, eigenstate));

            // Return the control qubit to the appropriate measurement basis.
            H(controlQubit);
        }

        adjoint auto
        controlled auto
        controlled adjoint auto
    }

    operation ContinuousPhaseEstimationIteration( oracle : ContinuousOracle, time : Double, theta : Double, eigenstate : Qubit[], controlQubit : Qubit)  : ()
    {
        body {
            if (time < 0.0) {
                // fail 'Oracle power cannot be negative.'
            }

            let inversionAngle = -(theta * time);

            // Prepare the control qubit.
            H(controlQubit);
            Rz(inversionAngle, controlQubit);

            (Controlled oracle)([controlQubit], (time, eigenstate));

            // Return the control qubit to the appropriate measurement basis.
            H(controlQubit);
        }

        adjoint auto
        controlled auto
        controlled adjoint auto
    }

    operation PrepAndMeasurePhaseEstImpl(wInv : Double, t : Double, op : ((Double, Double, Qubit) => ())) : Result {
        body {
            mutable datum = Zero;
            using (register = Qubit[1]) {
                op(t, wInv, register[0]);
                set datum = MResetZ(register[0]);
            }
            return datum;
        }
    }


}
