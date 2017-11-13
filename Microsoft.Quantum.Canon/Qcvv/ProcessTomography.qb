// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {
    open Microsoft.Quantum.Primitive;
    open Microsoft.Quantum.Extensions.Convert;
    open Microsoft.Quantum.Extensions.Math;

    operation MeasureAllZeroState(register : Qubit[]) : Result {
        body {
            // TODO: extract this into its own function.
            let nQubits = Length(register);
            mutable allZMeasurement = new Pauli[nQubits];
            for (idxQubit in 0..nQubits - 1) {
                set allZMeasurement[idxQubit] = PauliZ;
            }

            return Measure(allZMeasurement, register);
        }
    }

    operation MeasureIdentity(register : Qubit[]) : Result {
        body {
            return Zero;
        }
    }

    operation PrepareSingleQubitIdentity(qubit : Qubit) : () {
        body {
            ApplyPauli([RandomSingleQubitPauli()], [qubit]);
        }
    }

    operation PrepareIdentity(register : Qubit[]) : () {
        body {
            ApplyToEach(PrepareSingleQubitIdentity, register);
        }
    }

    operation EstimateFrequency(preparation : (Qubit[] => ()), measurement : (Qubit[] => Result), nQubits : Int, nMeasurements : Int) : Double {
        body {
            mutable nUp = 0;

            for (idxMeasurement in 0..nMeasurements - 1) {
                using (register = Qubit[nQubits]) {
                    preparation(register);
                    let result = measurement(register);
                    if (result == Zero) {
                        // NB!!!!! This reverses Zero and One to use conventions
                        //         common in the QCVV community. That is confusing
                        //         but is confusing with an actual purpose.
                        set nUp = nUp + 1;
                    }
                    // NB: We absolutely must reset here, since preparation()
                    //     and measurement() can each use randomness internally.
                    ApplyToEach(Reset, register);
                }
            }

            return ToDouble(nUp) / ToDouble(nMeasurements);
        }
    }

    /// summary:
    ///     Given two registers, prepares the maximally entangled state
    ///     |ß00><ß00| between each pair of qubits on the respective registers,
    ///     assuming that each register starts in the |0…0> state.
    operation PrepareEntangledState(left : Qubit[], right : Qubit[]) : () {
        body {
            if (Length(left) != Length(right)) {
                fail "Left and right registers must be the same length.";
            }

            for (idxQubit in 0..Length(left) - 1) {
                H(left[idxQubit]);
                (Controlled X)([left[idxQubit]], right[idxQubit]);
            }
        }

        adjoint auto
        controlled auto
        controlled adjoint auto
    }

    /// summary:
    ///     Prepares the Choi–Jamilkowski state for a given operation onto given reference
    ///     and target registers.
    operation PrepareChoiStateCA(op : (Qubit[] => () : Controlled, Adjoint), reference : Qubit[], target : Qubit[]) : () {
        body {
            PrepareEntangledState(reference, target);
            op(target);
        }

        adjoint auto
        controlled auto
        controlled adjoint auto
    }

}
