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

    // NOTE: the last qubit of 'reg' in this context is the auxillary qubit used in the Hadamard test.
    operation _ApplyLEOperationToRawRegister(op : (LittleEndian => Unit is Adj), target : Qubit[]) : Unit is Adj {
        op(LittleEndian(target));
    }

    operation _EstimateDerivativeWithParameterShift(
        inputEncoder : StateGenerator,
        gates : GateSequence,
        parameters : (Double[], Double[]),
        nQubits : Int,
        nMeasurements : Int
    ) : Double {
        return EstimateRealOverlapBetweenStates(
            _ApplyLEOperationToRawRegister(inputEncoder::Apply, _),
            _ApplyGates(Fst(parameters), gates, _),
            _ApplyGates(Snd(parameters), gates, _),
            nQubits, nMeasurements
        );
    }

    /// # Summary
    /// polymorphic classical/quantum gradient estimator
    ///
    /// # Input
    /// ## param
    /// circuit parameters
    ///
    /// ## gates
    /// sequence of gates in the circuits
    ///
    /// ## sg
    /// generates quantum encoding of a subject sample (either simulated or true)
    ///
    /// ## measCount
    /// number of true quantum measurements to estimate probabilities.
    /// IMPORTANT: measCount==0 implies simulator deployment
    ///
    /// # Output
    /// the gradient
    ///
    operation EstimateGradient(
        gates : GateSequence,
        param : Double[],
        sg : StateGenerator,
        nMeasurements : Int
    )
    : (Double[]) {
        //Synopsis: Suppose (param,gates) define Circ0
        //Suppose (param1,gates1) define Circ1 that implements one-gate derivative of Circ0
        //The expectation derivative is then 2 Re[<Circ1 psi|\Pi_1|Circ0 psi>] =
        // Re[<Circ1 psi|Id|Circ0 psi>] - Re[<Circ1 psi|Z \otimes Id|Circ0 psi>]
        //We observe SEE THEORY that for (Circ1)=(Circ0)' ,  Re[<Circ1 psi|Circ0 psi>]==0
        //Thus we are left to compute Re[<Circ1 psi|Z \otimes Id|Circ0 psi>] =
        // 1 - 1/2 < (Z \otimes Id) Circ0 psi - Circ1 psi | (Z \otimes Id) Circ0 psi - Circ1 psi>
        //i.e., 1 - HadamardTestResultHack(Circ1,[Z],Circ0)


        //Now, suppose a gate at which we differentiate is the (Controlled R(\theta))([k0,k1,...,kr],[target])
        //and we want a unitary description of its \theta-derivative. It can be written as
        // 1/2 {(Controlled R(\theta'))([k0,k1,...,kr],[target]) -  (Controlled Z)([k1,...,kr],[k0])(Controlled R(\theta'))([k0,k1,...,kr],[target])}
        mutable grad = ConstantArray(Length(param), 0.0);
        let nQubits = MaxI(NQubitsRequired(gates), sg::NQubits);

        for (gate in gates!) {
            let paramShift = (param + [0.0])
                // Shift the corresponding parameter.
                w/ gate::Index <- (param[gate::Index] + PI());

            // NB: This the *antiderivative* of the bracket
            let newDer = _EstimateDerivativeWithParameterShift(
                sg, gates, (param, paramShift), nQubits, nMeasurements
            );
            if (IsEmpty(gate::Span::ControlIndices)) {
                //uncontrolled gate
                set grad w/= gate::Index <- grad[gate::Index] + newDer;
            } else {
                //controlled gate
                let controlledShift = paramShift
                    w/ gate::Index <- (param[gate::Index] + 3.0 * PI());
                //Assumption: any rotation R has the property that R(\theta+2 Pi)=(-1).R(\theta)
                // NB: This the *antiderivative* of the bracket
                let newDer1 = _EstimateDerivativeWithParameterShift(
                    sg, gates, (param, controlledShift), nQubits, nMeasurements
                );
                set grad w/= gate::Index <- (grad[gate::Index] + 0.5 * (newDer - newDer1));
            }
        }
        return grad;

    }

}
