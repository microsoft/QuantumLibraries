namespace Microsoft.Quantum.MachineLearning {
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Canon;
	open Microsoft.Quantum.Convert;

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
				let circEnc = NoisyInputEncoder(effectiveTolerance, samp);
				set ret += [QubitProbPhysical(circEnc, param, gates, nQubits, nMeasurements)];
			}
		}

		return ret;
	} //EstimateClassificationProbabilitiesClassicalData

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
		return InferredLabels(probs, bias);
	}


}
