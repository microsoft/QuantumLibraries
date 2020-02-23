// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.MachineLearning {
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Preparation;
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
            if (idxNegative >= Length(coefficients)) {
                fail $"Cannot set the phase at index {idxNegative}, only {Length(coefficients)} coefficients were provided.";
            }
            let coefficient = coefficients[idxNegative];
            set ret w/= idxNegative <- ComplexPolar(coefficient::Magnitude, 0.0);
        }
        return ret;
    }

    function _NegativeLocations(cNegative: Int, coefficients : ComplexPolar[]) : Int[] {
        mutable negLocs = new Int[0];
        for ((idx, coefficient) in Enumerated(coefficients)) {
            if (AbsD(coefficient::Argument - PI()) < 1E-9) {
                set negLocs += [idx];
            }
        }
        return Length(negLocs) > cNegative ? negLocs[...cNegative - 1] | negLocs;
    }

    /// Do special processing on the first cNegative entries
    operation _ReflectAboutNegativeCoefficients(
        negLocs : Int[],
        coefficients : ComplexPolar[],
        reg: LittleEndian
    )
    : Unit is Adj + Ctl {
        for (idxNegative in negLocs) {
            ReflectAboutInteger(idxNegative, reg);
        }
    }

    /// # Summary
    /// Returns the number of qubits required to encode a particular feature
    /// vector.
    ///
    /// # Input
    /// ## sample
    /// A sample feature vector to be encoded into a qubit register.
    ///
    /// # Output
    /// The size required to encode `sample` into a qubit register, expressed
    /// as a number of qubits.
    function FeatureRegisterSize(sample : Double[]) : Int {
        return Ceiling(Lg(IntAsDouble(Length(sample))));
    }

    /// # Summary
    /// Given a set of coefficients and a tolerance, returns a state preparation
    /// operation that prepares each coefficient as the corresponding amplitude
    /// of a computational basis state, up to the given tolerance.
    ///
    /// # Input
    /// ## tolerance
    /// The approximation tolerance to be used in encoding the given
    /// coefficients into an input state.
    /// ## coefficients
    /// The coefficients to be encoded into an input state.
    ///
    /// # Output
    /// A state preparation operation that prepares the given coefficients
    /// as an input state on a given register.
    function ApproximateInputEncoder(tolerance : Double, coefficients : Double[])
    : StateGenerator {
        //First quantize the coefficients: for a coef x find such y*tolerance, where y is integer and |x-y*tolerance| \neq tolerance/2
        let nCoefficients = Length(coefficients);
        mutable complexCoefficients = new ComplexPolar[Length(coefficients)];
        mutable cNegative = 0;
        for ((idx, coef) in Enumerated(coefficients)) {
            mutable magnitude = coef;
            if (tolerance > 1E-9) {
                set magnitude = tolerance * IntAsDouble(Round(coefficients[idx] / tolerance)); //quantization
            }
            mutable ang = 0.0;
            if (magnitude < 0.0) {
                set cNegative += 1;
                set magnitude = -magnitude;
                set ang = PI();
            }
            set complexCoefficients w/= idx <- ComplexPolar(magnitude, ang);
        }

        // Check if we can apply the explicit two-qubit case.
        if (_CanApplyTwoQubitCase(coefficients)) {
            return StateGenerator(2, _ApplyTwoQubitCase(coefficients, _));
        }

        let nQubits = FeatureRegisterSize(coefficients);

        // If not, we may be able to use a special protocol in the case that
        // there are only a few negative coefficients.
        // Here, by a "few," we mean fewer than the number of qubits required
        // to encode features.
        if ((cNegative > 0) and (IntAsDouble(cNegative) < Lg(IntAsDouble(Length(coefficients))) + 1.0)) {
            let negLocs = _NegativeLocations(cNegative, complexCoefficients);
            return StateGenerator(
                nQubits,
                BoundCA([
                    // Prepare the state disregarding the sign of negative components.
                    _CompileApproximateArbitraryStatePreparation(
                        tolerance, _Unnegate(negLocs, complexCoefficients), nQubits
                    ),
                    // Reflect about the negative coefficients to apply the negative signs
                    // at the end.
                    _ReflectAboutNegativeCoefficients(negLocs, complexCoefficients, _)
                ])
            );
        }

        // Finally, we fall back to arbitrary state preparation.
        return StateGenerator(
            nQubits,
            _CompileApproximateArbitraryStatePreparation(tolerance, complexCoefficients, nQubits)
        );
    }

    /// # Summary
    /// Given a set of coefficients and a tolerance, returns a state preparation
    /// operation that prepares each coefficient as the corresponding amplitude
    /// of a computational basis state.
    ///
    /// # Input
    /// ## coefficients
    /// The coefficients to be encoded into an input state.
    ///
    /// # Output
    /// A state preparation operation that prepares the given coefficients
    /// as an input state on a given register.
    function InputEncoder(coefficients : Double[])
    : StateGenerator {
        //default implementation, does not respect sparcity
        mutable complexCoefficients = new ComplexPolar[Length(coefficients)];
        for ((idx, coefficient) in Enumerated(coefficients)) {
            set complexCoefficients w/= idx <- ComplexPolar(
                coefficient >= 0.0
                ? (coefficient, 0.0)
                | (-coefficient, PI())
            );
        }
        if (_CanApplyTwoQubitCase(coefficients)) {
            return StateGenerator(2, _ApplyTwoQubitCase(coefficients, _));
        }
        //this is preparing the state almost exactly so far
        return StateGenerator(
            FeatureRegisterSize(coefficients),
            ApproximatelyPrepareArbitraryState(1E-12, complexCoefficients, _)
        );
    }

}
