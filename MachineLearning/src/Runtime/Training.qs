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
    ) : SequentialModel {
        mutable bestSoFar = SequentialModel([-1E12], -2.0);
        mutable bestValidation = Length(samples) + 1;

        let features = Mapped(_Features, samples);
        let labels = Mapped(_Label, samples);

        for (idxStart in 0..(Length(parameterSource) - 1)) {
            Message($"Beginning training at start point #{idxStart}...");
            let ((h, m), proposedUpdate) = StochasticTrainingLoop(
                gates, SequentialModel(parameterSource[idxStart], 0.0),
                samples, options, trainingSchedule, 1
            );
            let probsValidation = EstimateClassificationProbabilitiesClassicalData(
                options::Tolerance, features, validationSchedule, nQubits,
                gates, proposedUpdate::Parameters, options::NMeasurements
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
                set bestSoFar = proposedUpdate;
            }

        }
        return bestSoFar;
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
        let options = DefaultTrainingOptions()
            w/ LearningRate <- learningRate
            w/ Tolerance <- tolerance
            w/ MinibatchSize <- miniBatchSize
            w/ NMeasurements <- nMeasurements
            w/ MaxEpochs <- maxEpochs;

        return (TrainSequentialClassifier(
            nQubits, gateSequence, parameterSource, samples,
            options, sch, schValidate
        ))!;
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
        miniBatch : LabeledSample[],
        options : TrainingOptions,
        param : Double[], gates : GateSequence
    )
    : (Double, Double[]) {
        mutable batchGradient = ConstantArray(Length(param), 0.0);

        for (samp in miniBatch) {
            mutable err = IntAsDouble(samp::Label);
            if (err < 1.0) {
                set err = -1.0; //class 0 misclassified to class 1; strive to reduce the probability
            }
            let grad = EstimateGradientFromClassicalSample(
                options::Tolerance, param, gates, samp::Features,
                options::NMeasurements
            );
            for (ip in 0..(Length(param) - 1)) {
                // GradientClassicalSample actually computes antigradient, but err*grad corrects it back to gradient
                set batchGradient w/= ip <- (batchGradient[ip] + options::LearningRate * err * grad[ip]);
            }

        }
        let updatedParameters = Mapped(PlusD, Zip(param, batchGradient));
        // TODO:REVIEW: Ok to interpret utility as size of the overall move?
        return (SquaredNorm(batchGradient), updatedParameters);
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
    operation OneStochasticTrainingEpoch(
        samples: LabeledSample[],
        schedule: SamplingSchedule, periodScore: Int,
        options : TrainingOptions,
        model : SequentialModel, gates: GateSequence,
        h0: Int, m0: Int
    )
    : ((Int, Int), SequentialModel) {
        let HARDCODEDunderage = 3; // 4/26 slack greater than 3 is not recommended


        mutable hBest = h0;
        mutable mBest = m0;
        mutable bestSoFar = model;

        let pls = ClassificationProbabilitiesClassicalData(
            samples, schedule, model::Parameters, gates, options::NMeasurements
        );
        let missLocations = MissLocations(schedule, pls, model::Bias);

        //An epoch is just an attempt to update the parameters by learning from misses based on LKG parameters
        let minibatches = Mapped(Subarray(_, samples), Chunks(options::MinibatchSize, missLocations));
        for (minibatch in minibatches) {
            let (utility, updatedParameters) = OneStochasticTrainingStep(
                minibatch, options, bestSoFar::Parameters, gates
            );
            if (utility > 0.0000001) {
                // There has been some good parameter update.
                // Check if it actually improves things, and if so,
                // commit it.
                let plsCurrent = ClassificationProbabilitiesClassicalData(samples, schedule, updatedParameters, gates, options::NMeasurements);
                let updatedBias = _UpdatedBias(
                    plsCurrent, model::Bias, options::Tolerance
                );
                let (h1, m1) = TallyHitsMisses(plsCurrent,  updatedBias);
                if (m1 < mBest + HARDCODEDunderage) {
                    //we allow limited non-greediness
                    if (m1 < mBest) {
                        set hBest = h1;
                        set mBest = m1;
                        set bestSoFar = SequentialModel(updatedParameters, updatedBias);
                    }
                }

            }

        }
        return ((hBest, mBest), bestSoFar);
    }

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
        gates: GateSequence,
        model : SequentialModel,
        samples: LabeledSample[],
        options : TrainingOptions,
        schedule: SamplingSchedule,
        periodScore: Int
    )
    : ((Int, Int), SequentialModel) {
        //const
        let relFuzz = 0.01;
        let pls = ClassificationProbabilitiesClassicalData(
            samples, schedule, model::Parameters, gates, options::NMeasurements
        );
        mutable bestSoFar = model
            w/ Bias <- _UpdatedBias(pls, model::Bias, options::Tolerance);
        mutable (hBest, mBest) = TallyHitsMisses(pls, model::Bias);
        mutable current = bestSoFar;

        //reintroducing learning rate heuristics
        mutable lrate = options::LearningRate;
        mutable batchSize = options::MinibatchSize;

        // Keep track of how many times a bias update has stalled out.
        mutable nStalls = 0;

        for (ep in 1..options::MaxEpochs) {
            let ((h1, m1), proposedUpdate) = OneStochasticTrainingEpoch(
                samples, schedule, periodScore,
                options
                    w/ LearningRate <- lrate
                    w/ MinibatchSize <- batchSize,
                current, gates, hBest, mBest
            );
            if (m1 < mBest) {
                set hBest = h1;
                set mBest = m1;
                set bestSoFar = proposedUpdate;
                if (IntAsDouble(mBest) / IntAsDouble(mBest + hBest) < options::Tolerance) { //Terminate based on tolerance
                    return ((hBest, mBest), bestSoFar);
                }
                set nStalls = 0; //Reset the counter of consequtive noops
                set lrate = options::LearningRate;
                set batchSize = options::MinibatchSize;
            }

            if (
                    NearlyEqualD(current::Bias, proposedUpdate::Bias) and _AllNearlyEqualD(current::Parameters, proposedUpdate::Parameters)
            ) {
                set nStalls += 1;
                // If we're more than halfway through our maximum allowed number of stalls,
                // exit early with the best we actually found.
                if (nStalls > options::MaxStalls) {
                    return ((hBest, mBest), bestSoFar); //Too many non-steps. Continuation makes no sense
                }

                // Otherwise, heat up the learning rate and batch size.
                set batchSize = nStalls; //batchSize + 1; //Try to fuzz things up with smaller batch count
                //and heat up  a bit
                set lrate *= 1.25;

                // If we stalled out, we'll also randomly rescale our parameters
                // and bias before updating.
                if (nStalls > options::MaxStalls / 2) {
                    set current = SequentialModel(
                        ForEach(_RandomlyRescale(relFuzz, _), proposedUpdate::Parameters),
                        _RandomlyRescale(relFuzz, proposedUpdate::Bias)
                    );
                }
            } else {
                // If we learned successfully this iteration, reset the number of
                // stalls so far.
                set nStalls = 0; //Reset the counter of consequtive noops
                set lrate = options::LearningRate;
                set batchSize = options::MinibatchSize;

                // Since we didn't stall out, we can set the parameters and bias
                // as normal, without randomizing.
                set current = proposedUpdate;
            }
        }

        return ((hBest, mBest), bestSoFar);
    }

}
