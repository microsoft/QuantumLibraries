// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {
    open Microsoft.Quantum.Primitive;

    /// <summary>
    ///     Cascade of multiply-controlled NOT gate with one or more control qubit and multiple target qubits. 
    /// </summary>
    /// <param name = "qs"> Quantum register of n qubits, to which an n qubit fanout gate is applied: qs[0] is control, all other qubits are targets. 
    operation MultiControlMultiTargetCNOT (controls : Qubit[], targets : Qubit[]) : () {
        body { 
            ApplyToEachAC(
                (Controlled X)(controls, _), targets
            );
        }
        adjoint self
        controlled auto
        controlled adjoint auto
    }

    /// <summary>
    /// Create a ket state, i.e., the state 1/sqrt(2)(|0...0> + |1..1>) on n qubits. 
    /// </summary>
    /// <param name = "qs"> Quantum register of n qubits, initially assumed to be in state |0...0> and on which the ket state is created. 
    operation CatState (qs : Qubit[]) : () {
        body { 
            if (Length(qs) < 2) { 
                fail "must have at least 2 qubits to create a cat state";
            }
            H(qs[0]);
            MultiControlMultiTargetCNOT([qs[0]], qs[1..Length(qs)-1]);
        }
    }

    // Functions to test the multiply controlled NOT gate and the cat state

    // FIXME: These can be turned into calls to ConstArray once generics support lands..
    function MakeAllXArray( length : Int ) : Pauli[] {
        mutable arr = new Pauli[length];
        for( i in 0 .. length - 1 ) { set arr[i] = PauliX; }
        return arr;
    }

    function MakeZZArray( length : Int, a : Int, b : Int ) : Pauli[] {
        mutable arr = new Pauli[length];
        for( i in 0 .. length - 1 ) { 
            if (i == a || i == b) { set arr[i] = PauliZ; } 
            else { set arr[i] = PauliI; }
        }
        return arr;
    }

    operation TestCatState (n : Int) : () {
        body {
            using( ancillas = Qubit[n] ) {
                CatState(ancillas);
                Assert(MakeAllXArray(n), ancillas, Zero, "Error: Cat state must be invariant under the X..X operator");
                for (i in 1..(n-1)) {
                    Assert(MakeZZArray(n, 0, i ), ancillas, Zero, "Error: Cat state must be invariant under all I..Z..Z..I operators");
                }
            }
        }
    }

    operation TestMultiControlMultiTargetCNOT (c : Int, n : Int) : () { 
        body {
            using ( controls = Qubit[c] ) {
                using ( targets = Qubit[n] ) {
                    for (i in 0..c-1) {
                        X(controls[i]);
                    }
                    MultiControlMultiTargetCNOT(controls, targets);
                    for (i in 0..(c-1)) {
                        AssertProb([PauliZ], [controls[i]], One, 1.0, "Error: Probability of measuring this qubit in One should be 1.0", 1e-10);
                    }
                    for (i in 0..(n-1)) {
                        AssertProb([PauliZ], [targets[i]], One, 1.0, "Error: Probability of measuring this qubit in One should be 1.0", 1e-10);
                    }
                }
            }
        }
    }
}
