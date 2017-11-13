// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {
    open Microsoft.Quantum.Primitive;

	/// summary:
	///     Teleportation transfers 1 qubit by encoding it into a 2-bit message,
	///     using an entangled pair of qubits.
	/// remarks:
    ///	    Always returns source qubit to |0>.
	///
	///     The circuit first creates an EPR pair between the target qubit and
	///     an ancilla qubit that gets allocated inside the function. Then a
	///     Bell measurement between the source qubit and one half of the EPR
	///     pair is performed. Finally, depending on the 4 possible outcomes of
	///     the Bell measurement, a correction is performed to restore the state
	///     in the target qubit.
	/// params:
	///     source: A single qubit representing the state to be teleported.
	///     target: A single qubit initially in the |0> state onto which
	///         given state is to be teleported.
	/// citation:
    ///     Nielsen & Chuang, CUP 2000, Section 1.3.6, http://doi.org/10.1017/CBO9780511976667
    operation Teleportation (source : Qubit, target : Qubit) : () {
        body {
            // Get a temporary qubit for the Bell pair.
            using (ancillaRegister = Qubit[1]) {
                let ancilla = ancillaRegister[0];
            
                // Create a Bell pair between the temporary qubit and the target.
                Assert([PauliZ], [target], Zero, "Error: target qubit must be initialized in zero state");
                H(ancilla);
                CNOT(ancilla, target);
                Assert([PauliZ; PauliZ], [ancilla; target], Zero, "Error: EPR state must be eigenstate of ZZ");
                Assert([PauliX; PauliX], [ancilla; target], Zero, "Error: EPR state must be eigenstate of XX");
            
                // Perform the Bell measurement and the correction necessary to
				// reconstruct the input state as the target state.
                CNOT(source, ancilla);
                H(source);
                AssertProb([PauliZ], [source], Zero, 0.5, "Error: All outcomes of the Bell measurement must be equally likely", 1e-5);
                AssertProb([PauliZ], [ancilla], Zero, 0.5, "Error: All outcomes of the Bell measurement must be equally likely", 1e-5);

                // Note that MResetZ makes sure that source is returned to zero state
				// so that we can deallocate it.
                if (M(source) == One)  { Z(target); }
				if (M(ancilla) == One) { X(target); }
            }
        }
    }


}
