// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {
    open Microsoft.Quantum.Primitive;

    /// <summary>
    /// Superdense coding transfers 2 classical bits by encoding them into 1 qubit, using 1 EPR pair ("2c=1q+1e")
    /// </summary>
    /// <param name = "a", "b"> Two Boolean values that represent the 2 input bits to the superdense coding protocol. </param>
    /// <remarks>  The circuit first creates an EPR pair between two ancilla qubits. Depending on the value of two classical 
    ///     bits then one out of 4 possible Bell states is created by applying a local transformation to just one half of the
    ///     EPR pair. Finally, a Bell measurement is applied to decode the two bits of classical information from the state.
    ///     [ Nielsen & Chuang, CUP 2000, Section 2.3, http://doi.org/10.1017/CBO9780511976667 ]
    /// </remarks>

    operation SuperdenseCoding (a : Bool, b : Bool) : (Result, Result) {
        body {
            mutable a_measured = Zero;
            mutable b_measured = Zero;
            // Get a temporary register for the Bell pair
            using (ancillas = Qubit[2]) {
                // Create an EPR state between the two ancilla qubits
                H(ancillas[0]);
                CNOT(ancillas[0], ancillas[1]);
                Assert([PauliZ; PauliZ], [ancillas[0]; ancillas[1]], Zero, "Error: EPR state must be eigenstate of ZZ");
                Assert([PauliX; PauliX], [ancillas[0]; ancillas[1]], Zero, "Error: EPR state must be eigenstate of XX");

                // Perform the superdense coding operation
                if a { Z(ancillas[0]); }
                if b { X(ancillas[0]); }

                // To decode, rotate into Bell basis and perform measurement
                CNOT(ancillas[0], ancillas[1]);
                H(ancillas[0]);

                // Perform checks if the measured result matches the inputs
                mutable a_expected = Zero;
                mutable b_expected = Zero;
                if a { set a_expected = One; }
                if b { set b_expected = One; }
                Assert([PauliZ], [ancillas[0]], a_expected, "Error: Superdense coding protocol must yield value of a in first qubit");
                Assert([PauliZ], [ancillas[1]], b_expected, "Error: Superdense coding protocol must yield value of b in second qubit");
                set a_measured = M(ancillas[0]);
                set b_measured = M(ancillas[1]);
            }
            // Return the result of the measurement
            return (a_measured, b_measured);
        }
    }

    operation SuperdenseCodingTestExFail () : () {
        body { 
            mutable res = (Zero, Zero);
            set res = SuperdenseCoding(true, true);
            set res = SuperdenseCoding(true, false);
            set res = SuperdenseCoding(false, true);
            set res = SuperdenseCoding(false, false);
        }
    }

}
