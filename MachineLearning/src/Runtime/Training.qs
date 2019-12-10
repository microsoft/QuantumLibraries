namespace Microsoft.Quantum.MachineLearning {
    open Microsoft.Quantum.Math;
    open Microsoft.Quantum.Logical;
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Convert;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Optimization;

    function _MisclassificationRate(probabilities : Double[], labels : Int[], bias : Double) : Double {
        let proposedLabels = InferredLabels(bias, probabilities);
        return IntAsDouble(NMismatches(proposedLabels, labels)) / IntAsDouble(Length(probabilities));
    }

    /// # Summary
    /// Returns a bias value that leads to near-minimum misclassification score.
    function _UpdatedBias(labeledProbabilities: (Double, Int)[], bias: Double, tolerance: Double) : Double {
        mutable min1 = 1.0;
        mutable max0 = 0.0;

        // Find the range of classification probabilities for each class.
        for ((probability, label) in labeledProbabilities) {
            if (label == 1) {
                if (min1 > probability) {
                    set min1 = probability;
                }
            } else {
                if (max0 < probability) {
                    set max0 = probability;
                }
            }
        }

        // Exit early if we can find a perfect classification.
        if (max0 <= min1) {
            return 0.5 * (1.0 - max0 - min1);
        }

        // If we can't find a perfect classification, minimize to find
        // the best feasible bias.
        let optimum = LocalUnivariateMinimum(
            _MisclassificationRate(Mapped(Fst<Double, Int>, labeledProbabilities), Mapped(Snd<Double, Int>, labeledProbabilities), _),
            (0.5 - max0, 0.5 - min1),
            tolerance
        );
        return optimum::Coordinate;
    }

    operation TrainSequentialClassifier(
        nQubits: Int,
        gates: GateSequence,
        parameterSource: Double[][],
        samples: LabeledSample[],
        options : TrainingOptions,
        trainingSchedule: SamplingSchedule,
        validationSchedule: SamplingSchedule
    ) : (Double[], Double) {
        mutable retParam = [-1E12];
        mutable retBias = -2.0; //Indicates non-informative start
        mutable bestValidation = Length(samples) + 1;

        let features = Mapped(_Features, samples);
        let labels = Mapped(_Label, samples);

        for (idxStart in 0..(Length(parameterSource) - 1)) {
            Message($"Beginning training at start point #{idxStart}...");
            let ((h, m), (b, parpar)) = StochasticTrainingLoop(
                samples, trainingSchedule, trainingSchedule, options, 1,
                parameterSource[idxStart], gates, 0.0
            );
            let probsValidation = EstimateClassificationProbabilitiesClassicalData(
                options::Tolerance, features, validationSchedule, nQubits,
                gates, parpar, options::NMeasurements
            );
            // Find the best bias for the new classification parameters.
            let localBias = _UpdatedBias(
                Zip(probsValidation, Sampled(validationSchedule, labels)),
                0.0,
                options::Tolerance
            );
            let localPL = InferredLabels(localBias, probsValidation);
            let localMisses = NMismatches(localPL, Sampled(validationSchedule, labels));
            if (bestValidation > localMisses) {
                set bestValidation = localMisses;
                set retParam = parpar;
                set retBias = localBias;
            }

        }
        return (retParam, retBias);
    }

    /// # Summary
    /// Using a flat description of a classification model, find a good local optimum
    /// for the model parameters and a related calssification bias
    ///
    /// # Input
    /// ## nQubits
    /// the number of qubits used for data encoding
    ///
    /// ## gates
    /// flat characterization of  circuit  structure. Each element is [parameterIndex, pauliCode, targetQubit\,sequence of control qubits\]
    ///
    /// ## parameterSource
    /// an array of parameter arrays, to be used as SGD starting points
    ///
    /// ## trainingSet
    /// the set of training samples
    ///
    /// ## trainingLabels
    /// the set of training labels
    ///
    /// ## trainingSchedule
    /// defines a subset of training data actually used in the training process
    ///
    /// ## validatioSchedule
    /// defines a subset of training data used for validation and computation of the *bias*
    ///
    /// ## learningRate
    /// initial learning rate for stochastic gradient descent
    ///
    /// ## tolerance
    /// sufficient absolute precision of parameter updates
    ///
    /// ## learningRate
    /// initial learning rate for stochastic gradient descent
    ///
    /// ## miniBatchSize
    /// maximum size of SGD mini batches
    ///
    /// ## maxEpochs
    /// limit to the number of training epochs
    ///
    /// ## nMeasurenets
    /// number of the measurement cycles to be used for estimation of each probability
    ///
    /// # Output
    /// (Array of optimal parameters, optimal validation *bias*)
    ///
    operation TrainQcccSequential(nQubits: Int, gates: Int[][], parameterSource: Double[][], trainingSet: Double[][], trainingLabels: Int[], trainingSchedule: Int[][], validationSchedule: Int[][],
                                    learningRate: Double, tolerance: Double, miniBatchSize: Int, maxEpochs: Int, nMeasurements: Int) : (Double[],Double) {
        let samples = unFlattenLabeledSamples(trainingSet,trainingLabels);
        let sch = unFlattenSchedule(trainingSchedule);
        let schValidate = unFlattenSchedule(validationSchedule);
        let gateSequence = unFlattenGateSequence(gates);
        let options = TrainingOptions(learningRate, tolerance, miniBatchSize, nMeasurements, maxEpochs);

        return TrainSequentialClassifier(
            nQubits, gateSequence, parameterSource, samples,
            options, sch, schValidate
        );
    } //TrainQcccSequential

    /// # Summary
    /// attempts a single parameter update in the direction of mini batch gradient
    ///
    /// # Input
    /// ## miniBatch
    /// container of labeled samples in the mini batch
    ///
    /// ## param
    /// circuit parameters
    ///
    /// ## gates
    /// sequence of gates in the circuits
    ///
    /// ## lrate
    /// the learning rate
    ///
    /// ## measCount
    /// number of true quantum measurements to estimate probabilities.
    ///
    /// # Output
    /// (utility, (new)parameters) pair
    ///
    operation OneStochasticTrainingStep(
        tolerance: Double, miniBatch: LabeledSample[], param: Double[], gates: GateSequence,
        lrate: Double, measCount: Int
    ) : (Double, Double[]) {
        mutable upParam = new Double[Length(param)];
        mutable batchGradient = ConstantArray(Length(param), 0.0);

        for (samp in miniBatch) {
            mutable err = IntAsDouble(samp::Label);
            if (err < 1.0) {
                set err = -1.0; //class 0 misclassified to class 1; strive to reduce the probability
            }
            let grad = EstimateGradientFromClassicalSample(tolerance, param, gates, samp::Features, measCount);
            for (ip in 0..(Length(param) - 1)) {
                // GradientClassicalSample actually computes antigradient, but err*grad corrects it back to gradient
                set batchGradient w/= ip <- (batchGradient[ip] + lrate * err * grad[ip]);
            }

        }
        for (ip in 0..(Length(param)-1)) {
            set upParam w/= ip <- (param[ip] + batchGradient[ip]);
        }
        return (SquaredNorm(batchGradient), upParam); //TODO:REVIEW: Ok to interpret utility as size of the overall move?
    }


    /// # Summary
    /// Perform one epoch of circuit training on a subset of data samples to a quantum simulator
    ///
    /// # Input
    /// ## samples
    /// a container of available data samples
    ///
    /// ## sched
    /// a schedule of the data subset for this training loop
    ///
    /// ## schedScore
    /// defines a (possibly different) data subset on which accuracy scoring is performed
    ///
    /// ## periodScore
    /// number of blind gradient steps between scoring points (performance tool, set to 1 for best accuracy)
    ///
    /// ## miniBatchSize
    /// number of samples in a gradient mini batch
    ///
    /// ## param
    /// initial parameter vector
    ///
    /// ## gates
    /// sequence of gates in the circuit
    ///
    /// ## bias
    /// reserved for future use; originally - initial prediction bias
    ///
    /// ## lrate
    /// learning rate
    ///
    /// ## measCount
    /// number of true quantum measurements to estimate probabilities.
    ///
    operation OneStochasticTrainingEpoch(samples: LabeledSample[], sched: SamplingSchedule, schedScore: SamplingSchedule, periodScore: Int,
                    miniBatchSize: Int, param: Double[], gates: GateSequence, bias: Double, lrate: Double, tolerance: Double, measCount: Int,
                    h0: Int, m0: Int): ((Int,Int),(Double,Double[]))
    {
        let HARDCODEDunderage = 3; //4/26 slack greater than 3 is not recommended


        mutable hBest = h0;
        mutable mBest = m0;
        mutable biasBest = bias;

        let pls = ClassificationProbabilitiesClassicalData(samples, schedScore, param, gates, measCount);
        let (h2,m2) = TallyHitsMisses(pls,biasBest);
        let missLocations = MissLocations(schedScore, pls, biasBest);

        mutable paramBest = param;
        mutable paramCurrent = paramBest;
        mutable biasCurrent = biasBest;

        //An epoch is just an attempt to update the parameters by learning from misses based on LKG parameters
        for (ixLoc in 0..miniBatchSize..(Length(missLocations) - 1)) {
            let miniBatch = ExtractMiniBatch(miniBatchSize, ixLoc, missLocations, samples);
            let (utility,upParam) = OneStochasticTrainingStep(tolerance, miniBatch, paramCurrent, gates, lrate, measCount);
            if (AbsD(utility) > 0.0000001) {
                //There had been some parameter update
                if (utility > 0.0) { //good parameter update
                    set paramCurrent = upParam;
                    let plsCurrent = ClassificationProbabilitiesClassicalData(samples, schedScore, paramCurrent, gates, measCount);
                    set biasCurrent = _UpdatedBias(plsCurrent, bias, tolerance);
                    let (h1,m1) = TallyHitsMisses(plsCurrent,biasCurrent);
                    if (m1 < mBest + HARDCODEDunderage) {
                        //we allow limited non-greediness
                        if (m1 < mBest) {
                            set hBest = h1;
                            set mBest = m1;
                            set paramBest = paramCurrent;
                            set biasBest = biasCurrent;
                        }
                    } else {
                        //otherwise we scrap the parameter update
                        set paramCurrent = paramBest;
                        set biasCurrent = biasBest;
                    }
                }

            }

        }
        return ((hBest, mBest), (biasBest, paramBest));
    }

    //Make some oblivious gradien descent steps without checking the prediction quality
    operation OneUncontrolledStochasticTrainingEpoch(samples: LabeledSample[], sched: SamplingSchedule, schedScore: SamplingSchedule, periodScore: Int,
                    miniBatchSize: Int, param: Double[], gates: GateSequence, bias: Double, lrate: Double, tolerance: Double,  measCount: Int): ((Int,Int),(Double,Double[]))
    {
        let pls = ClassificationProbabilitiesClassicalData(samples, schedScore, param, gates, measCount);
        mutable biasBest = _UpdatedBias(pls, bias, tolerance);
        let (h0,m0) = TallyHitsMisses(pls,biasBest); // ClassificationScoreSimulated(samples, schedScore, param, gates, bias); //Deprecated
        mutable hCur = h0;
        mutable mCur = m0;
        let missLocations = MissLocations(schedScore, pls, biasBest);

        mutable paramBest = param;
        mutable paramCurrent = paramBest;
        mutable biasCurrent = biasBest;

        //An epoch is just an attempt to update the parameters by learning from misses based on LKG parameters
        for (ixLoc in 0..miniBatchSize..(Length(missLocations) - 1)) {
            let miniBatch = ExtractMiniBatch(miniBatchSize,ixLoc,missLocations,samples);
            let (utility,upParam) = OneStochasticTrainingStep(tolerance, miniBatch, paramCurrent, gates, lrate, measCount);
            if (AbsD(utility) > 0.0000001) {
                //There had been some parameter update
                if (utility > 0.0) { //good parameter update
                    set paramCurrent = upParam;
                    let plsCurrent = ClassificationProbabilitiesClassicalData(samples, schedScore, paramCurrent, gates, measCount);
                    set biasCurrent = _UpdatedBias(plsCurrent, bias, tolerance);
                    let (h1,m1) = TallyHitsMisses(plsCurrent,biasCurrent);
                    set hCur = h1;
                    set mCur = m1;
                }

            }

        }
        return ((hCur, mCur),(biasCurrent,paramCurrent));
    } //OneUncontrolledStochasticTrainingEpoch

    /// # Summary
    /// Randomly rescales an input to either grow or shrink by a given factor.
    operation _RandomlyRescale(scale : Double, value : Double) : Double {
        return value * (
            1.0 + scale * (Random([0.5, 0.5]) > 0 ? 1.0 | -1.0)
        );
    }

    /// # Summary
    /// Run a full circuit training loop on a subset of data samples
    ///
    /// # Input
    /// ## samples
    /// a container of available data samples
    ///
    /// ## sched
    /// a schedule of the data subset for this training loop
    ///
    /// ## schedScore
    /// defines a (possibly different) data subset on which accuracy scoring is performed
    ///
    /// ## periodScore
    /// number of blind gradient steps between scoring points (performance tool, set to 1 for best accuracy)
    ///
    /// ## miniBatchSize
    /// number of samples in a gradient mini batch
    ///
    /// ## param
    /// initial parameter vector
    ///
    /// ## gates
    /// sequence of gates in the circuit
    ///
    /// ## bias
    /// reserved for future use; originally - initial prediction bias
    ///
    /// ## lrate
    /// learning rate
    ///
    /// ## maxEpochs
    /// maximum number of epochs in this loop
    ///
    /// ## tol
    /// tolerance: acceptable misprediction rate in training
    ///
    /// ## measCount
    /// number of true quantum measurements to estimate probabilities.
    /// IMPORTANT: measCount==0 implies simulator deployment
    ///
    /// # Output
    /// ((no.hits,no.misses),(opt.bias,opt.parameters))
    ///
    operation StochasticTrainingLoop(
        samples: LabeledSample[], sched: SamplingSchedule, schedScore: SamplingSchedule,
        options : TrainingOptions,
        periodScore: Int,
        param: Double[], gates: GateSequence,
        bias: Double
    )
    : ((Int, Int),(Double, Double[])) {
        //const
        let relFuzz = 0.01;
        let pls = ClassificationProbabilitiesClassicalData(samples, schedScore, param, gates, options::NMeasurements);
        mutable biasBest = _UpdatedBias(pls, bias, options::Tolerance);
        mutable (hBest, mBest) = TallyHitsMisses(pls, biasBest);
        mutable paramBest = param;
        mutable paramCurrent = param;
        mutable biasCurrent = biasBest;

        //reintroducing learning rate heuristics
        mutable lrate = options::LearningRate;
        mutable batchSize = options::MinibatchSize;

        // Keep track of how many times a bias update has stalled out.
        mutable nStalls = 0;

        for (ep in 1..options::MaxEpochs) {
            let ((h1, m1), (upBias, upParam)) = OneStochasticTrainingEpoch(
                samples, sched, schedScore, periodScore,
                batchSize, paramCurrent, gates, biasCurrent, lrate,
                options::Tolerance, options::NMeasurements, hBest, mBest
            );
            if (m1 < mBest) {
                set hBest = h1;
                set mBest = m1;
                set paramBest = upParam;
                set biasBest = upBias;
                if (IntAsDouble(mBest) / IntAsDouble(mBest + hBest) < options::Tolerance) { //Terminate based on tolerance
                    return ((hBest, mBest), (biasBest, paramBest));
                }
                set nStalls = 0; //Reset the counter of consequtive noops
                set lrate = options::LearningRate;
                set batchSize = options::MinibatchSize;
            }

            if (NearlyEqualD(biasCurrent, upBias) and _AllNearlyEqualD(paramCurrent, upParam)) {
                set nStalls += 1;
                // If we're more than halfway through our maximum allowed number of stalls,
                // exit early with the best we actually found.
                if (nStalls > options::MaxStalls) {
                    return ((hBest, mBest), (biasBest, paramBest)); //Too many non-steps. Continuation makes no sense
                }

                // Otherwise, heat up the learning rate and batch size.
                set batchSize = nStalls; //batchSize + 1; //Try to fuzz things up with smaller batch count
                //and heat up  a bit
                set lrate *= 1.25;

                // If we stalled out, we'll also randomly rescale our parameters
                // and bias before updating.
                if (nStalls > options::MaxStalls / 2) {
                    set biasCurrent = _RandomlyRescale(relFuzz, upBias);
                    set paramCurrent = ForEach(_RandomlyRescale(relFuzz, _), upParam);
                }
            } else {
                // If we learned successfully this iteration, reset the number of
                // stalls so far.
                set nStalls = 0; //Reset the counter of consequtive noops
                set lrate = options::LearningRate;
                set batchSize = options::MinibatchSize;

                // Since we didn't stall out, we can set the parameters and bias
                // as normal, without randomizing.
                set paramCurrent = upParam;
                set biasCurrent = upBias;
            }
        }

        return ((hBest, mBest), (biasBest, paramBest));
    }

}
