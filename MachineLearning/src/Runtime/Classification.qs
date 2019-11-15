namespace Microsoft.Quantum.MachineLearning {
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Canon;

    
	/// # Summary
	/// Using a flat description of a classification model, assign estimated probability of the top class label
	/// to each vector in the test set
	///
	/// # Input
	/// ## nQubits
	/// the number of qubits used for data encoding
	///
	/// ## gates
	/// flat characterization of  circuit  structure. Each element is [parameterIndex, pauliCode, targetQubit\,sequence of control qubits\]
	///
	/// ## parameters
	/// an array of circuit parameters
	///
	/// ## testSet
	/// the set of vectors to be labeled
	///
	/// ## nMeasurenets
	/// number of the measurement cycles to be used for estimation of each probability
	///
	/// # Output
	/// Array of estimated probabilities of top class label (for each sample in the test set)
	///
	operation EstimateClassificationProbabilities(tolerance: Double, nQubits: Int, gates: Int[][], parameters: Double[], testSet: Double[][], nMeasurements: Int) : Double[]
	{
		let segSched = [0..1..Length(testSet)-1];
		return EstimateClassificationProbabilitiesClassicalData(tolerance, testSet, SamplingSchedule(segSched), nQubits, unFlattenGateSequence(gates), parameters, nMeasurements);
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
	/// ## testSet
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
    operation DoClassification(tolerance: Double, nQubits: Int, gates: Int[][], parameters: Double[], bias: Double, testSet: Double[][], nMeasurements: Int) : Int[]
	{
		let probs = EstimateClassificationProbabilities(tolerance, nQubits,gates,parameters,testSet,nMeasurements);
		return InferredLabels(probs, bias);
	}


}
