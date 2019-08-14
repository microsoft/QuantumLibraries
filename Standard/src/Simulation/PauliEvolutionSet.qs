// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Simulation {
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Convert;

    // Convention for GeneratorIndex = ((Int[],Double[]), Int[])
    // We index single Paulis as 0 for I, 1 for X, 2 for Y, 3 for Z.
    // We index Pauli strings with arrays of integers e.g. a = [3,1,1,2] for ZXXY.
    // We assume the zeroth element of Double[] is the angle of rotation
    // We index the qubits that Pauli strings act on with arrays of integers e.g. q = [2,4,5,8] for Z_2 X_4 X_5, Y_8
    // An example of a Pauli string GeneratorIndex is thus ((a,b), q)

    // Consider the Hamiltonian H = 0.1 XI + 0.2 IX + 0.3 ZY
    // Its GeneratorTerms are (([1],b),[0]), 0.1),  (([1],b),[1]), 0.2),  (([3,2],b),[0,1]), 0.3).

    /// # Summary
    /// Represents a dynamical generator as a set of simulatable gates and an
    /// expansion in the Pauli basis.
    ///
    /// See [Dynamical Generator Modeling](/quantum/libraries/data-structures#dynamical-generator-modeling)
    /// for more details.
    ///
    /// # Input
    /// ## generatorIndex
    /// A generator index to be represented as unitary evolution in the Pauli
    /// basis.
    /// ## delta
    /// A multiplier on the duration of time-evolution by the term referenced
    /// in `generatorIndex`.
    /// ## qubits
    /// Register acted upon by time-evolution operator.
    operation PauliEvolutionImpl (generatorIndex : GeneratorIndex, delta : Double, qubits : Qubit[]) : Unit
    is Adj + Ctl {
        let (idxPaulis, idxDoubles) = generatorIndex::Data;
        let pauliString = IntArrayAsPauliArray(idxPaulis);
        let op = Exp(pauliString, delta * idxDoubles[0], _);
        (RestrictedToSubregisterCA(op, generatorIndex::Subsystems))(qubits);
    }

    /// # Summary
    /// Represents a dynamical generator as a set of simulatable gates and an
    /// expansion in the Pauli basis.
    ///
    /// # Input
    /// ## generatorIndex
    /// A generator index to be represented as unitary evolution in the Pauli
    /// basis.
    ///
    /// # Output
    /// An `EvolutionUnitary` representing time-evolution by the term
    /// referenced in `generatorIndex.
    function PauliEvolutionFunction(generatorIndex : GeneratorIndex) : EvolutionUnitary {
        return EvolutionUnitary(PauliEvolutionImpl(generatorIndex, _, _));
    }


    /// # Summary
    /// Represents a dynamical generator as a set of simulatable gates and an
    /// expansion in the Pauli basis.
    ///
    /// # Output
    /// An `EvolutionSet` that maps a `GeneratorIndex` for the Pauli basis to
    /// an `EvolutionUnitary.
    ///
    /// # Remarks
    /// This is obtained by partial application of
    /// <xref:microsoft.quantum.simulation.paulievolutionfunction>.
    function PauliEvolutionSet() : EvolutionSet {
        return EvolutionSet(PauliEvolutionFunction);
    }

}
