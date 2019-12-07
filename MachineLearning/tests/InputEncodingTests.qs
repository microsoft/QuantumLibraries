namespace Microsoft.Quantum.MachineLearning.Tests {
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Math;
    open Microsoft.Quantum.Arithmetic;
    open Microsoft.Quantum.Preparation;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Diagnostics;
    open Microsoft.Quantum.MachineLearning;

    operation _ApplyToBareRegister(op : (LittleEndian => Unit is Adj), register : Qubit[]) : Unit is Adj {
        op(LittleEndian(register));
    }

    @Test("QuantumSimulator")
    operation CheckInputEncoderWithPositiveInputs() : Unit {
        let coefficients = [0.1, 0.2, 0.3, 0.4];
        let encoder = InputEncoder(coefficients);
        AssertOperationsEqualReferenced(2,
            _ApplyToBareRegister(PrepareArbitraryState(Mapped(ComplexPolar(_, 0.0), coefficients), _), _),
            _ApplyToBareRegister(encoder, _)
        );
    }

}
