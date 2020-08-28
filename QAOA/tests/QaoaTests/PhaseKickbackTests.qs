// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Tests {
    open Microsoft.Quantum.Logical;
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Diagnostics;
    open Microsoft.Quantum.Math;
    open Microsoft.Quantum.Qaoa;

    @Test("QuantumSimulator")
    operation RunPhaseKickbackTest() : Unit {

        let phaseExponent = 0.5;

		using (register = Qubit()) {
            RunPhaseKickback(phaseExponent, [register]);
            AssertQubit(Zero, register);
            Reset(register);
        }
    }

     @Test("QuantumSimulator")
    operation RunPhaseKickbackOneControlQubitTest() : Unit {

        let phaseExponent = 0.5;
        
        let complexZero = Complex(0.685125,-0.174941);
        let complexOne = Complex(0.685125,0.174941);

        using (register = Qubit()) {
            H(register);
            RunPhaseKickback(phaseExponent, [register]);
            AssertQubitIsInStateWithinTolerance((complexZero, complexOne), register, 1E-05);
            Reset(register);
        }
    }
}
