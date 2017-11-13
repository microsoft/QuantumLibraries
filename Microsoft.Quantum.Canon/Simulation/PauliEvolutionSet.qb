// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {
    open Microsoft.Quantum.Primitive;


    // Convention for GeneratorIndex = ((Int[],Double[]), Int[])
    // We index single Paulis as 0 for I, 1 for X, 2 for Y, 3 for Z.
    // We index Pauli strings with arrays of integers e.g. a = [3;1;1;2] for ZXXY.
    // We assume the zeroth element of Double[] is the angle of rotation
    // We index the qubits that Pauli strings act on with arrays of integers e.g. q = [2;4;5;8] for Z_2 X_4 X_5, Y_8 
    // An example of a Pauli string GeneratorIndex is thus ((a,b), q)

    // Consider the Hamiltonian H = 0.1 XI + 0.2 IX + 0.3 ZY
    // Its GeneratorTerms are (([1],b),[0]), 0.1),  (([1],b),[1]), 0.2),  (([3;2],b),[0;1]), 0.3).

    function IntToPauli(idx : Int) : Pauli {
        let paulis = [PauliI; PauliX; PauliY; PauliZ];
        return paulis[idx];
    }

    /// # Summary
    /// Converts an array of integers to an array of single-qubit Pauli operators.
    ///
    /// # Input
    /// ## ints
    /// Array of integers in the range `0..3 - 1`  to be converted to Pauli
    /// operators.
    ///
    /// # Output
    /// An array `paulis` of Pauli operators the same length as `ints` such
    /// that `paulis[idxPauli]` is equal to the element of
    /// `[PauliI; PauliX; PauliY; PauliZ]` given by `ints[idxPauli]`.
    function IntsToPaulis(ints : Int[]) : Pauli[] {
        let nInts = Length(ints);
        mutable paulis = new Pauli[nInts];
        for (idxInt in 0..nInts - 1){
            set paulis[idxInt] = IntToPauli(ints[idxInt]);
        }

        return paulis;
    }

    /// summary:
    ///     Represents a dynamical generator as a set of simulatable gates and
    ///     an expansion in terms of that basis.
    ///     Las parameter for number of terms
    operation PauliEvolutionImpl(generatorIndex : GeneratorIndex, delta : Double, qubits: Qubit[]) : ()
    {
        body {
            let ((idxPaulis, idxDoubles), idxQubits) = generatorIndex;
            let pauliString = IntsToPaulis(idxPaulis);

            let op = Exp(pauliString, delta * idxDoubles[0], _);
            (RestrictToSubregisterCA(op, idxQubits))(qubits);
        }
        adjoint auto
        controlled auto
        controlled adjoint auto
    }
    function PauliEvolutionFunction(generatorIndex : GeneratorIndex) : EvolutionUnitary
    {
        return EvolutionUnitary(PauliEvolutionImpl(generatorIndex, _, _));
    }

    function PauliEvolutionSet() : EvolutionSet
    {
        return EvolutionSet(PauliEvolutionFunction);
    }

}
