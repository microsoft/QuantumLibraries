// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Preparation {
    open Microsoft.Quantum.Arithmetic;
    open Microsoft.Quantum.Math;

    @Deprecated("Microsoft.Quantum.Preparation.PrepareSingleQubitPositivePauliEigenstate")
    operation PrepareQubit(basis : Pauli, qubit : Qubit) : Unit {
        PrepareSingleQubitPositivePauliEigenstate(basis, qubit);
    }

    @Deprecated("Microsoft.Quantum.Preparation.PurifiedMixedState")
    function QuantumROM(targetError: Double, coefficients: Double[])
    : ((Int, (Int, Int)), Double, ((LittleEndian, Qubit[]) => Unit is Adj + Ctl)) {
        let preparation = PurifiedMixedState(targetError, coefficients);
        return (
            preparation::Requirements!,
            preparation::Norm,
            preparation::Prepare
        );
    }

    @Deprecated("Microsoft.Quantum.Preparation.PurifiedMixedStateRequirements")
    function QuantumROMQubitCount(targetError: Double, nCoeffs: Int)
    : (Int, (Int, Int)) {
        return (PurifiedMixedStateRequirements(targetError, nCoeffs))!;
    }

    @Deprecated("Microsoft.Quantum.Preparation.PrepareArbitraryStateCP")
    operation PrepareArbitraryState(coefficients : ComplexPolar[], qubits : LittleEndian)
    : Unit is Adj + Ctl {
        PrepareArbitraryStateCP(coefficients, qubits);
    }

    @Deprecated("Microsoft.Quantum.Preparation.PrepareArbitraryStateCP")
    function StatePreparationComplexCoefficients (coefficients : ComplexPolar[]) : (LittleEndian => Unit is Adj + Ctl) {
        return PrepareArbitraryStateCP(coefficients, _);
    }

    @Deprecated("Microsoft.Quantum.Preparation.PrepareArbitraryStateD")
    function StatePreparationPositiveCoefficients (coefficients : Double[])
    : (LittleEndian => Unit is Adj + Ctl) {
        let nCoefficients = Length(coefficients);
        mutable coefficientsComplexPolar = new ComplexPolar[nCoefficients];

        for (idx in 0 .. nCoefficients - 1) {
            set coefficientsComplexPolar w/= idx <- ComplexPolar(AbsD(coefficients[idx]), 0.0);
        }

        return PrepareArbitraryStateCP(coefficientsComplexPolar, _);
    }
}
