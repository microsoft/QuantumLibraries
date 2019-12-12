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

    /// WARNING: the downstream EstimateFrequencyA counts the frequency of Zero

    operation measureLastQubit(nQubits : Int): (Qubit[] => Result) {
        let paulis = ConstantArray(nQubits, PauliI) w/ (nQubits - 1) <- PauliZ;
        return Measure(paulis, _);
    }

    operation _endToEndPreparation(enc: (LittleEndian => Unit is Adj + Ctl), parameters: Double[], gates: GateSequence, reg: Qubit[]): Unit is Adj
    {
        enc(LittleEndian(reg));
        _ApplyGates(parameters, gates, reg);
    }

    operation endToEndPreparation(enc: (LittleEndian => Unit is Adj + Ctl), parameters: Double[], gates: GateSequence) : (Qubit[] => Unit is Adj)
    {
        return _endToEndPreparation(enc,parameters, gates, _);
    }

    function collectNegativeLocs(cNegative: Int, coefficients : ComplexPolar[]) : Int[]
    {
        mutable negLocs = ConstantArray(cNegative, -1);
        mutable nlx = 0;
        for (idx in 0 .. Length(coefficients) - 1)
        {
            let (r,a) = (coefficients[idx])!;
            if (AbsD(a - PI()) <  1E-9) {
                if (nlx < cNegative)
                {
                    set negLocs w/= nlx <- idx;
                    set nlx = nlx+1;
                }
            }
        }
        return negLocs;
    } //collectNegativeLocs

    // NOTE: the last qubit of 'reg' in this context is the auxillary qubit used in the Hadamard test.
    operation _endToEndHTcircuit(enc: (LittleEndian => Unit is Adj + Ctl), param1 : Double[], gates1: GateSequence, param2 : Double[], gates2: GateSequence, reg: Qubit[]): Unit is Adj + Ctl {
        let L = Length(reg) - 1;
        let g1 = _ApplyGates(param1,gates1,_);
        let g2 = _ApplyGates(param2,gates2,_);

        enc(LittleEndian(reg[0..(L-1)]));
        within {
            H(Tail(reg));
        } apply {
            (Controlled g1) ([reg[L]], reg[0..(L-1)]);
            within {
                X(Tail(reg));
            } apply {
                (Controlled g2) ([reg[L]], reg[0..(L-1)]);
                (Controlled Z)  ([reg[L]], reg[(L-1)]);
            }
        }
    }

    operation endToEndHTcircuit(enc: (LittleEndian => Unit is Adj + Ctl),param1 : Double[], gates1: GateSequence, param2 : Double[], gates2: GateSequence) : (Qubit[] => Unit is Adj) {
        return _endToEndHTcircuit(enc,param1, gates1, param2, gates2, _);
    }

    operation HardamardTestPhysical(enc2: (LittleEndian => Unit is Adj + Ctl), param1 : Double[], gates1: GateSequence, param2 : Double[], gates2: GateSequence, nQubits: Int, nMeasurements : Int): Double
    {
        return 1.0-EstimateFrequencyA(endToEndHTcircuit(enc2,param1,gates1,param2,gates2),measureLastQubit(nQubits), nQubits, nMeasurements);
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
    operation EstimateGradient(param : Double[], gates: GateSequence, sg: StateGenerator, nMeasurements : Int) : (Double[]) {
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
        let pC = Length(param);
        mutable grad = ConstantArray(pC, 0.0);
        mutable paramShift = param + [0.0];
        let nQubits = MaxI(NQubitsRequired(gates), sg::NQubits);

        for (gate in gates!) {
            set paramShift w/= gate::Index <- (param[gate::Index] + PI()); //Shift the corresponding parameter
            // NB: This the *antiderivative* of the bracket
            let newDer = 2.0 * HardamardTestPhysical(
                sg::Apply, param, gates, paramShift, gates, nQubits + 1, nMeasurements
            ) - 1.0;
            if (IsEmpty(gate::Span::ControlIndices)) {
                //uncontrolled gate
                set grad w/= gate::Index <- grad[gate::Index] + newDer;
            } else {
                //controlled gate
                set paramShift w/=gate::Index<-(param[gate::Index]+3.0 * PI());
                //Assumption: any rotation R has the property that R(\theta+2 Pi)=(-1).R(\theta)
                // NB: This the *antiderivative* of the bracket
                let newDer1 = 2.0 * HardamardTestPhysical(
                    sg::Apply, param, gates, paramShift, gates, nQubits + 1,
                    nMeasurements
                ) - 1.0;
                set grad w/= gate::Index <- (grad[gate::Index] + 0.5* (newDer - newDer1));
                set paramShift w/= gate::Index <-( param[gate::Index] + PI()); //unshift by 2 Pi (for debugging purposes)
            }
            set paramShift w/= gate::Index <- param[gate::Index]; //unshift this parameter
        }
        return grad;

    } //GradientHack


    /// # Summary
    /// computes stochastic gradient on one classical sample
    ///
    /// # Input
    /// ## param
    /// circuit parameters
    ///
    /// ## gates
    /// sequence of gates in the circuits
    ///
    /// ## sample
    /// sample vector as a raw array
    ///
    /// ## nMeasurements
    /// number of true quantum measurements to estimate probabilities
    ///
    /// # Output
    /// the gradient
    ///
    operation EstimateGradientFromClassicalSample(tolerance: Double, param : Double[], gates: GateSequence, sample: Double[], nMeasurements : Int) : (Double[]) {
        let nQubits = MaxI(FeatureRegisterSize(sample), NQubitsRequired(gates));
        let circEnc = NoisyInputEncoder(tolerance / IntAsDouble(Length(gates!)), sample);
        let sg = StateGenerator(nQubits, circEnc);
        return EstimateGradient(param, gates, sg, nMeasurements);
    }

    //Csharp-frendly adapter for gradient estimation
    //'gates' is a array of "flattened" controlled rotation defitions
    //each such definition is Int[no.controls+3] in the format [parameter index, Pauli index, target index <,control qubit indices>]
    //Pauli index is: 0 for I, 1 for X, 2 for y, 3 for Z
    //target index is the index of the target qubit of the rotation
    //Sequence of <control qubit indices> can be empty for uncontroled
    operation GradientClassicalSimulationAdapter(tolerance: Double, param : Double[], gates: Int[][], sample: Double[]) : (Double[])
    {

        return EstimateGradientFromClassicalSample(tolerance, param,unFlattenGateSequence(gates),sample,0);

    }

}
