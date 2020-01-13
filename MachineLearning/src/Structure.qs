// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.MachineLearning {
    open Microsoft.Quantum.Math;
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Canon;

    /// # Summary
    /// Returns the number of qubits required to apply a given sequential
    /// classifier.
    ///
    /// # Input
    /// ## structure
    /// The structure of a given sequential classifier.
    ///
    /// # Output
    /// The minimum size of a register on which the sequential classifier
    /// may be applied.
    function NQubitsRequired(structure : SequentialClassifierStructure)
    : Int {
        mutable nQubitsRequired = 0;
        for (gate in structure!) {
            set nQubitsRequired = Fold(
                MaxI, 0,
                gate::Span::ControlIndices + [
                    gate::Span::TargetIndex,
                    nQubitsRequired
                ]
            );
        }
        return nQubitsRequired;
    }

    /// # Summary
    /// Given the structure and parameterization of a sequential classifier,
    /// applies the classifier to a register of qubits.
    ///
    /// # Input
    /// ## structure
    /// Structure of the given sequential classifier.
    /// ## parameters
    /// A parameterization at which the given structure is applied.
    /// ## qubits
    /// A target register to which the classifier should be applied.
    operation ApplySequentialClassifier(
        structure : SequentialClassifierStructure,
        parameters : Double[],
        qubits : Qubit[]
    )
    : (Unit) is Adj + Ctl {
        for (gate in structure!) {
            if (gate::Index < Length(parameters)) {
                let input = (gate::Axis, parameters[gate::Index], qubits[gate::Span::TargetIndex]);
                if (IsEmpty(gate::Span::ControlIndices)) {
                    // Uncontrolled rotation of target
                    R(input);
                } else {
                    //TODO: should one validate the control indices first?
                    (Controlled R)(Subarray(gate::Span::ControlIndices, qubits), input);
                }
            }
        }
    }

    function _UncontrolledSpanSequence(idxsQubits : Int[]) : GateSpan[] {
        return Mapped(
            GateSpan(_, new Int[0]),
            idxsQubits
        );
    }

    function _CallFlipped<'TInput1, 'TInput2, 'TOutput>(
        fn : (('TInput1, 'TInput2) -> 'TOutput),
        y : 'TInput2, x : 'TInput1
    ) : 'TOutput {
        return fn(x, y);
    }

    function _Flipped<'TInput1, 'TInput2, 'TOutput>(
        fn : (('TInput1, 'TInput2) -> 'TOutput)
    ) : (('TInput2, 'TInput1) -> 'TOutput) {
        return _CallFlipped(fn, _, _);
    }

    function LocalRotationsLayer(nQubits : Int, axis : Pauli) : SequentialClassifierStructure {
        // [parameterIndex, pauliCode, targetQubit\,sequence of control qubits\]
        return SequentialClassifierStructure(Mapped(
            _Flipped(ControlledRotation(_, axis, _)),
            Enumerated(
                _UncontrolledSpanSequence(SequenceI(0, nQubits - 1))
            )
        ));
    }


    function PartialRotationsLayer(idxsQubits : Int[], axis : Pauli) : SequentialClassifierStructure {
        // [parameterIndex, pauliCode, targetQubit\,sequence of control qubits\]
        return SequentialClassifierStructure(Mapped(
            _Flipped(ControlledRotation(_, axis, _)),
            Enumerated(
                _UncontrolledSpanSequence(idxsQubits)
            )
        ));
    }

    function CyclicEntanglingLayer(nQubits : Int, axis : Pauli, stride : Int) : SequentialClassifierStructure {
        mutable rotations = new ControlledRotation[0];
        for (idxTarget in 0..nQubits - 1) {
            set rotations += [ControlledRotation(
                GateSpan(
                    idxTarget,
                    [(idxTarget + stride) % nQubits]
                ),
                axis, idxTarget
            )];
        }
        return SequentialClassifierStructure(rotations);
    }

    function CombinedStructure(layers : SequentialClassifierStructure[]) : SequentialClassifierStructure {
        mutable combined = (Head(layers))!;
        mutable offset = Length(combined);
        for (layer in Rest(layers)) {
            for (gate in layer!) {
                set combined += [gate w/ Index <- gate::Index + offset];
            }
            set offset += Length(layer!);
        }
        return SequentialClassifierStructure(combined);
    }

}