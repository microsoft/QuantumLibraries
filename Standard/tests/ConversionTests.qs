// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Tests {
    open Microsoft.Quantum.Diagnostics;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Convert;

    @Test("QuantumSimulator")
    function ResultArrayAsIntIsCorrect() : Unit {
        EqualityFactI(ResultArrayAsInt([Zero, Zero]), 0, $"Expected [Zero, Zero] to be represented by 0.");
        EqualityFactI(ResultArrayAsInt([One, Zero]), 1, $"Expected [One, Zero] to be represented by 1.");
        EqualityFactI(ResultArrayAsInt([Zero, One]), 2, $"Expected [Zero, One] to be represented by 2.");
        EqualityFactI(ResultArrayAsInt([One, One]), 3, $"Expected [One, One] to be represented by 3.");
    }

    @Test("QuantumSimulator")
    function BoolArrFromPositiveIntIsCorrect() : Unit {
        for number in 0 .. 100 {
            let bits = IntAsBoolArray(number, 9);
            let inte = BoolArrayAsInt(bits);
            EqualityFactI(inte, number, $"Integer converted to bit string and back should be identical");
        }

        let bits70 = [false, true, true, false, false, false, true, false];
        let number70 = BoolArrayAsInt(bits70);
        EqualityFactI(70, number70, $"Integer from 01000110 in little Endian should be 70");
    }

    @Test("QuantumSimulator")
    function IntAsBoolArrayIsCorrectForLargeInputs() : Unit {
        // Checks for regression of https://github.com/microsoft/QuantumLibraries/issues/503.
        let bits = IntAsBoolArray(1, 63);
        EqualityFactI(Length(bits), 63, $"IntAsBoolArray(1, 63) returned an array of the wrong length.");
        EqualityFactI(BoolArrayAsInt(bits), 1, $"IntAsBoolArray(1, 63) returned an array that did not represent 1: {bits}.");
    }

}
