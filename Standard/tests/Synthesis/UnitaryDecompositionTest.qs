// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Tests {
    open Microsoft.Quantum.Arithmetic;
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Convert;
    open Microsoft.Quantum.Diagnostics;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Math;
    open Microsoft.Quantum.Measurement;
    open Microsoft.Quantum.Synthesis;
    open Microsoft.Quantum.Random;

 

    @Test("QuantumSimulator")
    operation MyFirstQuantumTest () : Unit {
        AssertOperationsEqualInPlace(1, ApplyUnitary, X);
    }
}
