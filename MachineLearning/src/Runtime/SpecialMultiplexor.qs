// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {
	open Microsoft.Quantum.Arithmetic;
	open Microsoft.Quantum.Intrinsic;
	open Microsoft.Quantum.Arrays;
	open Microsoft.Quantum.Math;
	
	/// # Summary
	/// Applies a Pauli rotation conditioned on an array of qubits.
	///
	/// This applies the multiply-controlled unitary operation $U$ that performs
	/// rotations by angle $\theta_j$ about single-qubit Pauli operator $P$
	/// when controlled by the $n$-qubit number state $\ket{j}$.
	///
	/// $U = \sum^{2^n-1}_{j=0}\ket{j}\bra{j}\otimes e^{i P \theta_j}$.
	///
	/// # Input
	/// ## tolerance
	/// Coefficients under this tolerance level should be ignored
	/// ## coefficients
	/// Array of up to $2^n$ coefficients $\theta_j$. The $j$th coefficient
	/// indexes the number state $\ket{j}$ encoded in little-endian format.
	///
	/// ## pauli
	/// Pauli operator $P$ that determines axis of rotation.
	///
	/// ## control
	/// $n$-qubit control register that encodes number states $\ket{j}$ in
	/// little-endian format.
	///
	/// ## target
	/// Single qubit register that is rotated by $e^{i P \theta_j}$.
	///
	/// # Remarks
	/// `coefficients` will be padded with elements $\theta_j = 0.0$ if
	/// fewer than $2^n$ are specified.
	operation NoisyMultiplexPauli (tolerance: Double,coefficients : Double[], pauli : Pauli, control : LittleEndian, target : Qubit) : Unit
	{
		body (...)
		{
			if (pauli == PauliZ)
			{
				let op = NoisyMultiplexZ(tolerance, coefficients, control, _);
				op(target);
			}
			elif (pauli == PauliX)
			{
				let op = NoisyMultiplexPauli(tolerance,coefficients, PauliZ, control, _);
				ApplyWithCA(H, op, target);
			}
			elif (pauli == PauliY)
			{
				let op = NoisyMultiplexPauli(tolerance,coefficients, PauliX, control, _);
				ApplyWithCA(Adjoint S, op, target);
			}
			elif (pauli == PauliI)
			{
				NoisyApplyDiagonalUnitary(tolerance,coefficients, control);
			}
			else
			{
				fail $"MultiplexPauli failed. Invalid pauli {pauli}.";
			}
		}
	
		adjoint invert;
		controlled distribute;
		controlled adjoint distribute;
	}
	

	function significantReal(tol: Double, rg:Double[]):Bool
	{
		for(j in 0..(Length(rg)-1))
		{
			if (AbsD(rg[j])>tol)
			{
				return true;
			}
		}
		return false;
	}
	
	/// # Summary
	/// Applies a Pauli Z rotation conditioned on an array of qubits.
	/// 
	/// This applies the multiply-controlled unitary operation $U$ that performs
	/// rotations by angle $\theta_j$ about single-qubit Pauli operator $Z$
	/// when controlled by the $n$-qubit number state $\ket{j}$.
	///
	/// $U = \sum^{2^n-1}_{j=0}\ket{j}\bra{j}\otimes e^{i Z \theta_j}$.
	///
	/// # Input
	/// ## coefficients
	/// Array of up to $2^n$ coefficients $\theta_j$. The $j$th coefficient
	/// indexes the number state $\ket{j}$ encoded in little-endian format.
	///
	/// ## control
	/// $n$-qubit control register that encodes number states $\ket{j}$ in
	/// little-endian format.
	///
	/// ## target
	/// Single qubit register that is rotated by $e^{i P \theta_j}$.
	///
	/// # Remarks
	/// `coefficients` will be padded with elements $\theta_j = 0.0$ if
	/// fewer than $2^n$ are specified.
	///
	/// # References
	/// - Synthesis of Quantum Logic Circuits
	/// Vivek V. Shende, Stephen S. Bullock, Igor L. Markov
	/// https://arxiv.org/abs/quant-ph/0406176
	operation NoisyMultiplexZ (tolerance: Double, coefficients : Double[], control : LittleEndian, target : Qubit) : Unit
	{
		body (...)
		{
			// pad coefficients length at tail to a power of 2.
			let coefficientsPadded = Padded(-2 ^ Length(control!), 0.0, coefficients);
	
			if (Length(coefficientsPadded) == 1)
			{
				// Termination case
				if (AbsD(coefficientsPadded[0])> tolerance)
				{
					Exp([PauliZ], coefficientsPadded[0], [target]);
				}
			}
			else
			{
				// Compute new coefficients.
				let (coefficients0, coefficients1) = specialMultiplexZComputeCoefficients_(coefficientsPadded);
				NoisyMultiplexZ(tolerance,coefficients0, LittleEndian((control!)[0 .. Length(control!) - 2]), target);
				if (significantReal(tolerance,coefficients1))
				{
					CNOT((control!)[Length(control!) - 1], target);
					NoisyMultiplexZ(tolerance,coefficients1, LittleEndian((control!)[0 .. Length(control!) - 2]), target);
					CNOT((control!)[Length(control!) - 1], target);
				}
			}
		}
	
		adjoint invert;
	
		controlled (controlRegister, ...)
		{
			// pad coefficients length to a power of 2.
			let coefficientsPadded = Padded(2 ^ (Length(control!) + 1), 0.0, Padded(-2 ^ Length(control!), 0.0, coefficients));
			let (coefficients0, coefficients1) = specialMultiplexZComputeCoefficients_(coefficientsPadded);
			NoisyMultiplexZ(tolerance,coefficients0, control, target);
			if (significantReal(tolerance,coefficients1))
			{
				Controlled X(controlRegister, target);
				NoisyMultiplexZ(tolerance,coefficients1, control, target);
				Controlled X(controlRegister, target);
			}
		}
	
		controlled adjoint invert;
	}
	
	
	/// # Summary
	/// Applies an array of complex phases to numeric basis states of a register of qubits.
	/// 
	/// That is, this implements the diagonal unitary operation $U$ that applies a complex phase
	/// $e^{i \theta_j}$ on the $n$-qubit number state $\ket{j}$.
	///
	/// $U = \sum^{2^n-1}_{j=0}e^{i\theta_j}\ket{j}\bra{j}$.
	///
	/// TODO: REIMPLEMENT THIS along the Welch et Bocharov lines
	/// # Input
	/// ## tolerance
	/// Coefficients under this tolerance level should be ignored
	/// ## coefficients
	/// Array of up to $2^n$ coefficients $\theta_j$. The $j$th coefficient
	/// indexes the number state $\ket{j}$ encoded in little-endian format.
	///
	/// ## control
	/// $n$-qubit control register that encodes number states $\ket{j}$ in
	/// little-endian format.
	///
	/// # Remarks
	/// `coefficients` will be padded with elements $\theta_j = 0.0$ if
	/// fewer than $2^n$ are specified.
	///
	/// # References
	/// - Synthesis of Quantum Logic Circuits
	/// Vivek V. Shende, Stephen S. Bullock, Igor L. Markov
	/// https://arxiv.org/abs/quant-ph/0406176
	operation NoisyApplyDiagonalUnitary (tolerance: Double, coefficients : Double[], qubits : LittleEndian) : Unit
	{
		body (...)
		{
			if (IsEmpty(qubits!)) {
				fail "operation ApplyDiagonalUnitary -- Number of qubits must be greater than 0.";
			}
	
			// pad coefficients length at tail to a power of 2.
			let coefficientsPadded = Padded(-2 ^ Length(qubits!), 0.0, coefficients);
	
			// Compute new coefficients.
			let (coefficients0, coefficients1) = specialMultiplexZComputeCoefficients_(coefficientsPadded);
			NoisyMultiplexZ(tolerance,coefficients1, LittleEndian((qubits!)[0 .. Length(qubits!) - 2]), (qubits!)[Length(qubits!) - 1]);
	
			if (Length(coefficientsPadded) == 2)
			{
				
				// Termination case
				if (AbsD(coefficients0[0])>tolerance)
				{
					Exp([PauliI], 1.0 * coefficients0[0], qubits!);
				}
			}
			else
			{
				NoisyApplyDiagonalUnitary(tolerance,coefficients0, LittleEndian((qubits!)[0 .. Length(qubits!) - 2]));
			}
		}
	
		adjoint invert;
		controlled distribute;
		controlled adjoint distribute;
	}
	
	
	/// # Summary
	/// Implementation step of multiply-controlled Z rotations.
	/// # See Also
	/// - Microsoft.Quantum.Canon.MultiplexZ
	function specialMultiplexZComputeCoefficients_ (coefficients : Double[]) : (Double[], Double[])
	{
		let newCoefficientsLength = Length(coefficients) / 2;
		mutable coefficients0 = new Double[newCoefficientsLength];
		mutable coefficients1 = new Double[newCoefficientsLength];
	
		for (idxCoeff in 0 .. newCoefficientsLength - 1)
		{
			set coefficients0 w/= idxCoeff <- 0.5 * (coefficients[idxCoeff] + coefficients[idxCoeff + newCoefficientsLength]);
			set coefficients1 w/= idxCoeff <- 0.5 * (coefficients[idxCoeff] - coefficients[idxCoeff + newCoefficientsLength]);
		}
	
		return (coefficients0, coefficients1);
	}
	

	
	
}
	

