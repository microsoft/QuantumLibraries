namespace Microsoft.Quantum.MachineLearning {
    open Microsoft.Quantum.Math;
    open Microsoft.Quantum.Characterization;
    open Microsoft.Quantum.Arithmetic;
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Canon;
	open Microsoft.Quantum.Convert;

	operation EstimateClassificationProbabilityFromEncodedSample(
		encodedSample : StateGenerator,
		parameters: Double[],
		gates: GateSequence, nMeasurements : Int
	)
	: Double {
		return 1.0 - EstimateFrequencyA(
			endToEndPreparation(encodedSample::Apply, parameters,gates),
			measureLastQubit(encodedSample::NQubits),
			encodedSample::NQubits,
			nMeasurements
		);
	}

	operation EstimateClassificationProbabilityFromSample(tolerance: Double, parameters : Double[], gates: GateSequence, sample: Double[], nMeasurements: Int)
	: Double {
		let nQubits = FeatureRegisterSize(sample);
		let circEnc = NoisyInputEncoder(tolerance / IntAsDouble(Length(gates!)), sample);
		return EstimateClassificationProbabilityFromEncodedSample(
			StateGenerator(nQubits, circEnc), parameters, gates, nMeasurements
		);

	}

	/// # Summary
	/// Given a of classification probability and a bias, returns the
	/// label inferred from that probability.
	///
	/// # Input
	/// ## bias
	/// The bias between two classes, typically the result of training a
	/// classifier.
	/// ## probability
	/// A classification probabilities for a particular sample, typicaly
	/// resulting from estimating its classification frequency.
	///
	/// # Output
	/// The label inferred from the given classification probability.
	function InferredLabel(bias : Double, probability : Double) : Int {
		return probability + bias > 0.5 ? 1 | 0;
	}

	/// # Summary
	/// Given an array of classification probabilities and a bias, returns the
	/// label inferred from each probability.
	///
	/// # Input
	/// ## bias
	/// The bias between two classes, typically the result of training a
	/// classifier.
	/// ## probabilities
	/// An array of classification probabilities for a set of samples, typicaly
	/// resulting from estimating classification frequencies.
	///
	/// # Output
	/// The label inferred from each classification probability.
	function InferredLabels(bias : Double, probabilities : Double[]): Int[] {
		return Mapped(InferredLabel(bias, _), probabilities);
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
	/// ## nQubits
	/// number of qubits in the classification circuit
	///
	/// ## gates
	/// the sequence of gates in the circuit
	///
	/// ## param
	/// parameters of the circuits
	///
	/// ## measCount
	///
	/// # Output
	/// array of corresponding estimated probabilities of the top class label
	///
	operation EstimateClassificationProbabilitiesClassicalData(
		tolerance : Double, samples : Double[][], sched : SamplingSchedule,
		nQubits : Int, gates : GateSequence, param : Double[],
		nMeasurements : Int
	) : Double[] {
		let effectiveTolerance = tolerance / IntAsDouble(Length(gates!));
		mutable ret = new Double[0];
		for (rg in sched!) {
			for (ix in rg) {
				let samp = samples[ix];
				set ret += [EstimateClassificationProbabilityFromEncodedSample(
					StateGenerator(nQubits, NoisyInputEncoder(effectiveTolerance, samp)),
					param, gates, nMeasurements
				)];
			}
		}

		return ret;
	}

	/// # Summary
	/// Using a flat description of a classification model, assign estimated probability of top class label
	/// to each vector in the test set
	///
	/// # Input
	/// ## nQubits
	/// the number of qubits used for data encoding
	///
	/// ## gates
	/// Flattened representation of classifier structure. Each element is
	/// [parameterIndex, pauliCode, targetQubit, sequence of control qubits]
	///
	/// ## parameters
	/// an array of circuit parameters
	///
	/// ## samples
	/// the set of vectors to be labeled
	///
	/// ## bias
	/// top class bias
	///
	/// ## nMeasurenets
	/// number of the measurement cycles to be used for estimation of each probability
	///
	/// # Output
	/// Array of predicted class labels for each sample of the test set
	///
    operation DoClassification(tolerance: Double, nQubits: Int, gates: Int[][], parameters: Double[], bias: Double, samples : Double[][], nMeasurements: Int) : Int[] {
		let schedule = SamplingSchedule([0..Length(samples) - 1]);
		let sequence = unFlattenGateSequence(gates);
		let probs = EstimateClassificationProbabilitiesClassicalData(
			tolerance, samples, schedule, nQubits, sequence, parameters, nMeasurements
		);
		return InferredLabels(bias, probs);
	}

}
