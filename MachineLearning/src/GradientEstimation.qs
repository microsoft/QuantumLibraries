// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.MachineLearning {
    open Microsoft.Quantum.Math;
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Arithmetic;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Convert;
    open Microsoft.Quantum.Diagnostics;
    open Microsoft.Quantum.Preparation;
    open Microsoft.Quantum.Characterization;

    // NOTE: the last qubit of 'reg' in this context is the auxiliary qubit used in the Hadamard test.
    operation _ApplyLEOperationToRawRegister(op : (LittleEndian => Unit is Adj), target : Qubit[]) : Unit is Adj {
        op(LittleEndian(target));
    }

    operation _EstimateDerivativeWithParameterShift(
        inputEncoder : StateGenerator,
        model : SequentialModel,
        parameters : (Double[], Double[]),
        nQubits : Int,
        nMeasurements : Int
    ) : Double {
        return EstimateRealOverlapBetweenStates(
            _ApplyLEOperationToRawRegister(inputEncoder::Prepare, _),
            ApplySequentialClassifier(model w/ Parameters <- Fst(parameters), _),
            ApplySequentialClassifier(model w/ Parameters <- Snd(parameters), _),
            nQubits, nMeasurements
        );
    }

    /// # Summary
    /// Estimates the training gradient for a sequential classifier at a
    /// particular model and for a given encoded input.
    ///
    /// # Input
    /// ## model
    /// The sequential model whose gradient is to be estimated.
    /// ## encodedInput
    /// An input to the sequential classifier, encoded into a state preparation
    /// operation.
    /// ## nMeasurements
    /// The number of measurements to use in estimating the gradient.
    ///
    /// # Output
    /// An estimate of the training gradient at the given input and model
    /// parameters.
    ///
    /// # Remarks
    /// This operation uses a Hadamard test and the parameter shift technique
    /// together to estimate the gradient.
    operation EstimateGradient(
        model : SequentialModel,
        encodedInput : StateGenerator,
        nMeasurements : Int
    )
    : (Double[]) {
        // Synopsis: Suppose (param,gates) define Circ0
        // Suppose (param1,gates1) define Circ1 that implements one-gate derivative of Circ0
        // The expectation derivative is then 2 Re[<Circ1 psi|\Pi_1|Circ0 psi>] =
        //  Re[<Circ1 psi|Id|Circ0 psi>] - Re[<Circ1 psi|Z \otimes Id|Circ0 psi>]
        // We observe SEE THEORY that for (Circ1)=(Circ0)' ,  Re[<Circ1 psi|Circ0 psi>]==0
        // Thus we are left to compute Re[<Circ1 psi|Z \otimes Id|Circ0 psi>] =
        //  1 - 1/2 < (Z \otimes Id) Circ0 psi - Circ1 psi | (Z \otimes Id) Circ0 psi - Circ1 psi>
        // i.e., 1 - HadamardTestResultHack(Circ1,[Z],Circ0)


        // Now, suppose a gate at which we differentiate is the (Controlled R(\theta))([k0,k1,...,kr],[target])
        // and we want a unitary description of its \theta-derivative. It can be written as
        // 1/2 {(Controlled R(\theta'))([k0,k1,...,kr],[target]) -  (Controlled Z)([k1,...,kr],[k0])(Controlled R(\theta'))([k0,k1,...,kr],[target])}
        mutable grad = ConstantArray(Length(model::Parameters), 0.0);
        let nQubits = MaxI(NQubitsRequired(model), encodedInput::NQubits);

        for (gate in model::Structure) {
            let paramShift = (model::Parameters + [0.0])
                // Shift the corresponding parameter.
                w/ gate::ParameterIndex <- (model::Parameters[gate::ParameterIndex] + PI());

            // NB: This the *antiderivative* of the bracket
            let newDer = _EstimateDerivativeWithParameterShift(
                encodedInput, model, (model::Parameters, paramShift), nQubits, nMeasurements
            );
            if (IsEmpty(gate::ControlIndices)) {
                //uncontrolled gate
                set grad w/= gate::ParameterIndex <- grad[gate::ParameterIndex] + newDer;
            } else {
                //controlled gate
                let controlledShift = paramShift
                    w/ gate::ParameterIndex <- (model::Parameters[gate::ParameterIndex] + 3.0 * PI());
                // Assumption: any rotation R has the property that R(\theta + 2 Pi) = (-1) R(\theta).
                // NB: This the *antiderivative* of the bracket
                let newDer1 = _EstimateDerivativeWithParameterShift(
                    encodedInput, model, (model::Parameters, controlledShift), nQubits, nMeasurements
                );
                set grad w/= gate::ParameterIndex <- (grad[gate::ParameterIndex] + 0.5 * (newDer - newDer1));
            }
        }
        return grad;

    }

}
