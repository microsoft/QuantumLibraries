// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {
    open Microsoft.Quantum.Primitive;

	/// <summary>
    ///     Given two operations, asserts that they act identically for all input states.
    ///     This assertion is implemented by using the Choi–Jamiłkowski isomorphism to reduce
    ///     the assertion to one of a qubit state assertion on two entangled registers.
    ///     Thus, this operation needs only a single call to each operation being tested,
    ///     but requires twice as many qubits to be allocated.
    ///     This assertion can be used to ensure, for instance, that an optimized version of an
    ///     operation acts identically to its naïve implementation, or that an operation
    ///     which acts on a range of non-quantum inputs agrees with known cases.
    /// </summary>
    /// <remarks>
    ///     This operation requires that the operation modeling the expected behavior is
    ///     adjointable, so that the inverse can be performed on the target register alone.
    ///     Formally, one can specify a transpose operation, which relaxes this requirement,
    ///     but the transpose operation is not in general physically realizable for arbitrary
    ///     quantum operations and thus is not included here as an option.
    /// </remarks>
    /// <param name="actual">Operation to be tested.</param>
    /// <param name="expected">Operation defining the expected behavior for the operation under test.</param>
    /// <param name="nQubits">Number of qubits to pass to each operation.</param>
    operation AssertOperationsEqualReferenced(actual : (Qubit[] => ()), expected : (Qubit[] => () : Adjoint), nQubits : Int) : () {
        body {
            // Prepare a reference register entangled with the target register.
            using (reference = Qubit[nQubits]) {
                using (target = Qubit[nQubits]) {
                    // FIXME: change to use PrepareChoiState/With/Bind, but that depends on generics since the arguments take
                    //        (Qubit[], Qubit[]), not Qubit[].
                    PrepareEntangledState(reference, target);
                    actual(target);
                    (Adjoint expected)(target);
                    (Adjoint PrepareEntangledState)(reference, target);

					let ket0 = (Complex(1.0, 0.0), Complex(0.0, 0.0));

                    ApplyToEach(AssertQubitState(ket0, _, 1e-5), reference);
                    ApplyToEach(AssertQubitState(ket0, _, 1e-5), target);
                }
                ResetAll(reference);
            }
        }
    }

}
