// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Diagnostics {
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Convert;
    open Microsoft.Quantum.Math;

    /// # Summary
    /// Asserts that the phase of an equal superposition state has the expected value.
    ///
    /// # Description
    /// This operation asserts that the phase $\phi$ of a quantum state
    /// that may be expressed as
    /// $\frac{e^{i t}}{\sqrt{2}}(e^{i\phi}\ket{0} + e^{-i\phi}\ket{1})$
    /// for some arbitrary real $t$ has the expected value.
    ///
    /// # Input
    /// ## expected
    /// The expected value of $\phi \in (-\pi,\pi]$.
    ///
    /// ## qubit
    /// The qubit that stores the expected state.
    ///
    /// ## tolerance
    /// Absolute tolerance on the difference between actual and expected.
    ///
    /// # Example
    /// The following assert succeeds:
    /// `qubit` is in state $\ket{\psi}=e^{i 0.5}\sqrt{1/2}\ket{0}+e^{i 0.5}\sqrt{1/2}\ket{1}$;
    /// - `AssertPhase(0.0, qubit, 10e-10);`
    ///
    /// `qubit` is in state $\ket{\psi}=e^{i 0.5}\sqrt{1/2}\ket{0}+e^{-i 0.5}\sqrt{1/2}\ket{1}$;
    /// - `AssertPhase(0.5, qubit, 10e-10);`
    ///
    /// `qubit` is in state $\ket{\psi}=e^{-i 2.2}\sqrt{1/2}\ket{0}+e^{i 0.2}\sqrt{1/2}\ket{1}$;
    /// - `AssertPhase(-1.2, qubit, 10e-10);`
    operation AssertPhase(expected : Double, qubit : Qubit, tolerance : Double) : Unit {
        let exptectedProbX = Cos(expected) * Cos(expected);
        let exptectedProbY = Sin(-1.0 * expected + PI() / 4.0) * Sin(-1.0 * expected + PI() / 4.0);
        AssertMeasurementProbability([PauliZ], [qubit], Zero, 0.5, $"AssertPhase failed. Was not given a uniform superposition.", tolerance);
        AssertMeasurementProbability([PauliY], [qubit], Zero, exptectedProbY, $"AssertPhase failed. PauliY Zero basis did not give probability {exptectedProbY}.", tolerance);
        AssertMeasurementProbability([PauliX], [qubit], Zero, exptectedProbX, $"AssertPhase failed. PauliX Zero basis did not give probability {exptectedProbX}.", tolerance);
    }

}
