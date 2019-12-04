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
				let op = ApproximatelyMultiplexZ(tolerance, coefficients, control, _);
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
				ApproximatelyApplyDiagonalUnitary(tolerance, coefficients, control);
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
	

