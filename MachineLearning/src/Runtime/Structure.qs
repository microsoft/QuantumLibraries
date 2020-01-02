// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.MachineLearning {

    open Microsoft.Quantum.Arrays;

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

    function LocalRotationsLayer(nQubits : Int, axis : Pauli) : GateSequence {
        // [parameterIndex, pauliCode, targetQubit\,sequence of control qubits\]
        return GateSequence(Mapped(
            _Flipped(ControlledRotation(_, axis, _)),
            Enumerated(
                _UncontrolledSpanSequence(SequenceI(0, nQubits - 1))
            )
        ));
    }


    function PartialRotationsLayer(idxsQubits : Int[], axis : Pauli) : GateSequence {
        // [parameterIndex, pauliCode, targetQubit\,sequence of control qubits\]
        return GateSequence(Mapped(
            _Flipped(ControlledRotation(_, axis, _)),
            Enumerated(
                _UncontrolledSpanSequence(idxsQubits)
            )
        ));
    }

    function CyclicEntanglingLayer(nQubits : Int, axis : Pauli, stride : Int) : GateSequence {
        mutable rotations = new ControlledRotation[0];
        for (idxTarget in 0..nQubits - 1) {
            set rotations += [ControlledRotation(
                GateSpan(
                    idxTarget,
                    [idxTarget + stride % nQubits]
                ),
                axis, idxTarget
            )];
        }
        return GateSequence(rotations);
    }

    function CombinedGateSequence(layers : GateSequence[]) : GateSequence {
        mutable combined = (Head(layers))!;
        mutable offset = Length(combined);
        for (layer in Rest(layers)) {
            for (gate in layer!) {
                set combined += [gate w/ Index <- gate::Index + offset];
            }
            set offset += Length(layer!);
        }
        return GateSequence(combined);
    }

}