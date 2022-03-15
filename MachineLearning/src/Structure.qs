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
    /// ## model
    /// The model for a given sequential classifier.
    ///
    /// # Output
    /// The minimum size of a register on which the sequential classifier
    /// may be applied.
    function NQubitsRequired(model : SequentialModel)
    : Int {
        mutable lastQubitIndex = -1; // Need to return 0 if there are no gates
        for gate in model::Structure {
            set lastQubitIndex = Fold(
                MaxI, lastQubitIndex,
                gate::ControlIndices + [gate::TargetIndex]
            );
        }
        return lastQubitIndex + 1;
    }

    /// # Summary
    /// Given the structure and parameterization of a sequential classifier,
    /// applies the classifier to a register of qubits.
    ///
    /// # Input
    /// ## model
    /// The sequential model describing the parameters and structure of the
    /// classifier to be applied.
    /// ## qubits
    /// A target register to which the classifier should be applied.
    operation ApplySequentialClassifier(
        model : SequentialModel,
        qubits : Qubit[]
    )
    : (Unit) is Adj + Ctl {
        for gate in model::Structure {
            if gate::ParameterIndex < Length(model::Parameters) {
                let input = (gate::Axis, model::Parameters[gate::ParameterIndex], qubits[gate::TargetIndex]);
                if IsEmpty(gate::ControlIndices) {
                    // Uncontrolled rotation of target
                    R(input);
                } else {
                    Controlled R(Subarray(gate::ControlIndices, qubits), input);
                }
            }
        }
    }

    function _UncontrolledSpanSequence(idxsQubits : Int[]) : (Int, Int[])[] {
        return Zipped(
            idxsQubits,
            ConstantArray(Length(idxsQubits), [])
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

    /// # Summary
    /// Returns an array of uncontrolled (single-qubit) rotations along a given
    /// axis, with one rotation for each qubit in a register, parameterized by
    /// distinct model parameters.
    ///
    /// # Input
    /// ## nQubits
    /// The number of qubits acted on by the given layer.
    /// ## axis
    /// The rotation axis for each rotation in the given layer.
    ///
    /// # Output
    /// An array of controlled rotations about the given axis, one on each of
    /// `nQubits` qubits.
    function LocalRotationsLayer(nQubits : Int, axis : Pauli) : ControlledRotation[] {
        // [parameterIndex, pauliCode, targetQubit\,sequence of control qubits\]
        return Mapped(
            _Flipped(ControlledRotation(_, axis, _)),
            Enumerated(
                _UncontrolledSpanSequence(SequenceI(0, nQubits - 1))
            )
        );
    }


    /// # Summary
    /// Returns an array of single-qubit rotations along a given
    /// axis, parameterized by distinct model parameters.
    ///
    /// # Input
    /// ## idxsQubits
    /// Indices for the qubits to be used as the targets for each rotation.
    /// ## axis
    /// The rotation axis for each rotation in the given layer.
    ///
    /// # Output
    /// An array of controlled rotations about the given axis, one on each of
    /// `nQubits` qubits.
    function PartialRotationsLayer(idxsQubits : Int[], axis : Pauli) : ControlledRotation[] {
        return Mapped(
            _Flipped(ControlledRotation(_, axis, _)),
            Enumerated(
                _UncontrolledSpanSequence(idxsQubits)
            )
        );
    }

    /// # Summary
    /// Returns an array of singly controlled rotations along a given axis,
    /// arranged cyclically across a register of qubits, and parameterized by
    /// distinct model parameters.
    ///
    /// # Input
    /// ## nQubits
    /// The number of qubits acted on by the given layer.
    /// ## axis
    /// The rotation axis for each rotation in the given layer.
    /// ## stride
    /// The separation between the target and control indices for each rotation.
    ///
    /// # Output
    /// An array of two-qubit controlled rotations laid out cyclically across
    /// a register of `nQubits` qubits.
    ///
    /// # Example
    /// The following are equivalent:
    /// ```qsharp
    /// let layer = CyclicEntanglingLayer(3, PauliX, 2);
    /// let layer = [
    ///     ControlledRotation((0, [2]), PauliX, 0),
    ///     ControlledRotation((1, [0]), PauliX, 1),
    ///     ControlledRotation((2, [1]), PauliX, 2)
    /// ];
    /// ```
    function CyclicEntanglingLayer(nQubits : Int, axis : Pauli, stride : Int) : ControlledRotation[] {
        mutable rotations = [];
        for idxTarget in 0..nQubits - 1 {
            set rotations += [ControlledRotation(
                (
                    idxTarget,
                    [(idxTarget + stride) % nQubits]
                ),
                axis, idxTarget
            )];
        }
        return rotations;
    }

    /// # Summary
    /// Given one or more layers of controlled rotations, returns a single
    /// layer with model parameter index shifted such that distinct layers
    /// are parameterized by distinct model parameters.
    ///
    /// # Input
    /// ## layers
    /// The layers to be combined.
    ///
    /// # Output
    /// A single layer of controlled rotations, representing the concatenation
    /// of all other layers.
    ///
    /// # Example
    /// The following are equivalent:
    /// ```qsharp
    /// let structure = CombinedStructure([
    ///     LocalRotationLayer(2, PauliY),
    ///     CyclicEntanglingLayer(3, PauliX, 2)
    /// ]);
    /// let structure = [
    ///     ControlledRotation((0, new Int[0]), PauliY, 0),
    ///     ControlledRotation((1, new Int[0]), PauliY, 1),
    ///     ControlledRotation((0, [2]), PauliX, 2),
    ///     ControlledRotation((1, [0]), PauliX, 3),
    ///     ControlledRotation((2, [1]), PauliX, 4)
    /// ];
    /// ```
    function CombinedStructure(layers : ControlledRotation[][]) : ControlledRotation[] {
        mutable combined = Head(layers);
        mutable offset = Length(combined);
        for layer in Rest(layers) {
            for gate in layer {
                set combined += [gate w/ ParameterIndex <- gate::ParameterIndex + offset];
            }
            set offset += Length(layer);
        }
        return combined;
    }

}
