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

    /// # Summary
    /// Get a list of all the classification probabilities. In the from of (prob1,label) pairs. THIS operation is IN DEPRECATION
    ///
    /// # Input
    /// ## samples
    /// a container of labeled samples
    ///
    /// ## sched
    /// a schedule to define a subset of samples
    ///
    /// ## param
    /// parameters of the circuits
    ///
    /// ## gates
    /// the sequence of gates in the circuit
    ///
    /// ## nMeasurements
    /// the maximum number of quantum measurements used in the probability estimation
    ///
    /// # Output
    /// TODO
    operation ClassificationProbabilitiesClassicalData(samples: LabeledSample[], sched: SamplingSchedule, param: Double[], gates: GateSequence, nMeasurements: Int):
        (Double,Int)[] {
        mutable N = IsEmpty(samples)
                    ? NQubitsRequired(gates)
                    | MaxI(NQubitsRequired(gates), FeatureRegisterSize(_Features(Head(samples))));
        mutable ret = new (Double, Int)[0];
        for (rg in sched!) {
            for (ix in rg) {
                let sample = samples[ix];
                //agnostic w.r.t. simulator (may still be simulable)
                let prob1 = EstimateClassificationProbabilityFromSample(1E-12, param, gates, sample::Features, nMeasurements);
                set ret += [(prob1, sample::Label)];
            }
        }

        return ret;
    }

    operation EstimateClassificationProbabilitiesClassicalDataAdapter(tolerance: Double, samples: Double[][], schedule: Int[][], nQubits: Int,  gates: Int[][], param: Double[], measCount: Int): Double[]
    {
        return EstimateClassificationProbabilitiesClassicalData(tolerance, samples, unFlattenSchedule(schedule), nQubits, unFlattenGateSequence(gates), param, measCount);
    }


    /// # Summary
    /// generate a flat list of sample indices where mispredictions occur
    ///
    /// # Input
    /// ## sched
    /// a sampling schedule
    ///
    /// ## pls
    /// a list of estimated probabilities with the corresponding class labels
    ///
    /// ## bias
    /// bias on record
    ///
    /// # Output
    /// the list of indices where mispredictions occur
    ///
    function MissLocations(sched : SamplingSchedule, pls : (Double, Int)[], bias: Double) : Int[] {
        mutable ret = new Int[0];
        mutable ir = 0;

        for (rg in sched!) {
            for (ix in rg) {
                let (prob1, lab) = pls[ir];
                set ir += 1;
                if (prob1 + bias > 0.5) {
                    if (lab < 1) {
                        set ret += [ix];
                    }
                } else {
                    if (lab > 0) {
                        set ret += [ix];
                    }
                }
            }
        }
        return ret;
    }

    /// # Summary
    /// C#-friendly adapter to misclassification tally
    ///
    /// # Input
    /// ## vectors
    /// data vectors in flat encoding
    ///
    /// ## labels
    /// array of corresponding class lables
    ///
    /// ## schedule
    /// flat representation of index subset on which the circuit is scored
    ///
    /// ## param
    /// circuit parameters
    ///
    /// ## gateStructure
    /// gate structure in flat representation
    ///
    /// ## bias
    /// prediction bias to be tested
    ///
    /// ## measCount
    /// maximum number of quantum measurements per estimation (measCount==0 implies simulator deployment)
    ///
    /// # Output
    /// the number of misclassifications
    ///
    operation MisclassificationScoreAdapter(vectors: Double[][], labels: Int[], schedule: Int[][], param: Double[], gateStructure: Int[][], bias: Double, measCount: Int) : Int {
        mutable misses = 0;
        let samples = unFlattenLabeledSamples(vectors,labels);
        let gates = unFlattenGateSequence(gateStructure);
        let sched = unFlattenSchedule(schedule);

        let pls = ClassificationProbabilitiesClassicalData(samples,sched,param,gates,measCount);
        let biasCurrent = _UpdatedBias(pls, bias, 0.01);
        let (h1,m1) = TallyHitsMisses(pls,biasCurrent);
        return m1;
    }

    /// # Summary
    /// Extract a mini batch of samples and wrap the batch as a LabeledSampleContainer
    ///
    /// # Input
    /// ## size
    /// desired number of samples in the mini batch
    ///
    /// ## ixLoc
    /// starting index for the batch in the list of locations
    ///
    /// ## locations
    /// list of indices of samples of interest
    ///
    /// ## samples
    /// the container to extract the samples from
    ///
    /// # Output
    /// the mini batched wrapped as a LabeledSampleContainer
    ///
    /// # Remarks
    /// the resulting mini batch can be occasionally shorter than the requested 'size'
    /// (when it falls on the tail end of the list of 'locations')
    ///
    @Deprecated("Microsoft.Quantum.Arrays.Chunks")
    function ExtractMiniBatch(size: Int, ixLoc: Int, locations: Int[], samples: LabeledSample[]): LabeledSample[] {
        let cnt = MinI(size, Length(locations) - ixLoc);
        mutable rgSamples = new LabeledSample[0];
        for (location in locations[ixLoc..ixLoc + cnt]) {
            set rgSamples += [samples[location]];
        }
        return rgSamples;
    }

}
