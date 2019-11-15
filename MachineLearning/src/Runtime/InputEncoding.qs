// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.MachineLearning {
    open Microsoft.Quantum.Convert;
    open Microsoft.Quantum.Math;
    open Microsoft.Quantum.Arithmetic;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Canon;

	function _CanApplyTwoQubitCase(datum: Double[]) : Bool {
		return((Length(datum)==4) and (Microsoft.Quantum.Math.AbsD(datum[0]*datum[3]-datum[1]*datum[2])< 1E-12) and (Microsoft.Quantum.Math.AbsD(datum[0])> 1E-4));
	}

	operation _ApplyTwoQubitCase(datum: Double[], reg: LittleEndian) : Unit is Adj + Ctl {
		let x = datum[1]/datum[0];
		let y = datum[2]/datum[0];
		// we now encoding [1,x,y,x*y]
		let ax = 2.0 * ArcTan(x);
		let ay = 2.0 * ArcTan(y);
		R(PauliY, ay, (reg!)[1]);
		R(PauliY, ax, (reg!)[0]);
	}

	function _Unnegate(negLocs: Int[], coefficients : ComplexPolar[]) : ComplexPolar[] {
		mutable ret = coefficients;
		for (idxNegative in negLocs) {
			let jx = negLocs[idxNegative];
			let coefficient = coefficients[jx];
			set ret w/= jx <- ComplexPolar(coefficient::Magnitude, 0.0);
		}
		return ret;
	}

	/// Do special processing on the first cNegative entries
	operation _EncodeSparseNegativeInput(cNegative: Int, tolerance: Double,coefficients : ComplexPolar[], reg: LittleEndian): Unit is Adj + Ctl
	{
		let negLocs = collectNegativeLocs(cNegative, coefficients);
		// Prepare the state disregarding the sign of negative components.
		NoisyPrepareArbitraryState(tolerance, _Unnegate(negLocs, coefficients), reg);
		// Reflect about the negative coefficients to apply the negative signs
		// at the end.
		for (ineg in 0..(cNegative - 1)) {
			let jx = negLocs[ineg];
			if (jx > -1) {
				ReflectAboutInteger(jx, reg); //TODO:REVIEW: this assumes that 2^Length(reg) is the minimal pad to Length(coefficients)
			}
		}
	}

	function NoisyInputEncoder(tolerance: Double,coefficients : Double[]) : (LittleEndian => Unit is Adj + Ctl) {
		//First quantize the coefficients: for a coef x find such y*tolerance, where y is integer and |x-y*tolerance| \neq tolerance/2
		let nCoefficients = Length(coefficients);
        mutable coefficientsComplexPolar = new ComplexPolar[nCoefficients];
        mutable cNegative = 0;
        for (idx in 0 .. nCoefficients - 1) {
			mutable coef = coefficients[idx];
			if (tolerance > 1E-9) {
				set coef = tolerance * IntAsDouble(Round(coefficients[idx] / tolerance)); //quantization
			}
			mutable ang = 0.0;
			if (coef < 0.0) {
				set cNegative += 1;
				set coef = -coef;
				set ang = PI();
			}
            set coefficientsComplexPolar w/= idx <- ComplexPolar(coef, ang);
        }

		// Check if we can apply the explicit two-qubit case.
        if (_CanApplyTwoQubitCase(coefficients)) {
			return _ApplyTwoQubitCase(coefficients, _);
		}
		// If not, we may be able to use a special protocol in the case that
		// there are only a few negative coefficients.
		// Here, by a "few," we mean fewer than the number of qubits required
		// to encode features.
		if ((cNegative > 0) and (IntAsDouble(cNegative) < Lg(IntAsDouble(Length(coefficients))) + 1.0)) {
			return _EncodeSparseNegativeInput(cNegative, tolerance, coefficientsComplexPolar, _); //TODO:MORE:ACCEPTANCE ("Wines" passing soi far)
		}
		
		// Finally, we fall back to arbitrary state preparation.
		return NoisyPrepareArbitraryState(tolerance, coefficientsComplexPolar, _);
	} //EncodeNoisyInput

	//TODO:REVIEW: Design consideration! The implicit qubit count must be read off from the state encoder, NOT from the gate sequence!

	/// Create amplitude encoding of an array of real-valued coefficients
	/// The vector of 'coefficients' does not have to be unitary
	function InputEncoder(coefficients : Double[]): (LittleEndian => Unit is Adj + Ctl) {
		//default implementation, does not respect sparcity
		let nCoefficients = Length(coefficients);
        mutable coefficientsComplexPolar = new ComplexPolar[nCoefficients];
        mutable allPositive = true;
        for (idx in 0 .. nCoefficients - 1) {
			mutable coef = coefficients[idx];
			mutable ang = 0.0;
			if (coef < 0.0)
			{
				set allPositive = false;
				set coef =  -coef;
				set ang =Microsoft.Quantum.Math.PI();
			}
            set coefficientsComplexPolar w/= idx<-ComplexPolar(coef,ang);
        }
        if (_CanApplyTwoQubitCase(coefficients)) {
			return _ApplyTwoQubitCase(coefficients,_);
		}
		return NoisyPrepareArbitraryState(1E-12, coefficientsComplexPolar, _); //this is preparing the state almost exactly so far
	}

}