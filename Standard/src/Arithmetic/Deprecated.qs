// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {
    open Microsoft.Quantum.Arithmetic;
    open Microsoft.Quantum.Measurement;
    open Microsoft.Quantum.Warnings;

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Arithmetic.ApplyReversedOpLE".
    operation ApplyReversedOpLittleEndian(op : (LittleEndian => Unit), register : BigEndian) : Unit {
        _Renamed("ApplyReversedOpLittleEndian", "ApplyReversedOpLE");
        ApplyReversedOpLE(op, register);
    }

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Arithmetic.ApplyReversedOpLEA".
    operation ApplyReversedOpLittleEndianA(op : (LittleEndian => Unit is Adj), register : BigEndian) : Unit {
        body (...) {
            _Renamed("ApplyReversedOpLittleEndianA", "ApplyReversedOpLEA");
            ApplyReversedOpLEA(op, register);
        }
        adjoint auto;
    }

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Arithmetic.ApplyReversedOpLEC".
    operation ApplyReversedOpLittleEndianC(op : (LittleEndian => Unit is Ctl), register : BigEndian) : Unit {
        body (...) {
            _Renamed("ApplyReversedOpLittleEndianC", "ApplyReversedOpLEC");
            ApplyReversedOpLEC(op, register);
        }
        controlled auto;
    }

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Arithmetic.ApplyReversedOpLEC".
    operation ApplyReversedOpLittleEndianCA(op : (LittleEndian => Unit is Adj + Ctl), register : BigEndian) : Unit {
        body (...) {
            _Renamed("ApplyReversedOpLittleEndianCA", "ApplyReversedOpLECA");
            ApplyReversedOpLECA(op, register);
        }
        adjoint auto;
        controlled auto;
        controlled adjoint auto;
    }

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Arithmetic.ApplyReversedOpBE".
    operation ApplyReversedOpBigEndian(op : (BigEndian => Unit), register : LittleEndian) : Unit {
        _Renamed("ApplyReversedOpBigEndian", "ApplyReversedOpBE");
        ApplyReversedOpBE(op, register);
    }

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Arithmetic.ApplyReversedOpBEA".
    operation ApplyReversedOpBigEndianA(op : (BigEndian => Unit is Adj), register : LittleEndian) : Unit {
        body (...) {
            _Renamed("ApplyReversedOpBigEndianA", "ApplyReversedOpBEA");
            ApplyReversedOpBEA(op, register);
        }
        adjoint auto;
    }
    
    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Arithmetic.ApplyReversedOpBEC".
    operation ApplyReversedOpBigEndianC(op : (BigEndian => Unit is Ctl), register : LittleEndian) : Unit {
        body (...) {
            _Renamed("ApplyReversedOpBigEndianC", "ApplyReversedOpBEC");
            ApplyReversedOpBEC(op, register);
        }
        controlled auto;
    }

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Arithmetic.ApplyReversedOpBEC".
    operation ApplyReversedOpBigEndianCA(op : (BigEndian => Unit is Adj + Ctl), register : LittleEndian) : Unit {
        body (...) {
            _Renamed("ApplyReversedOpBigEndianCA", "ApplyReversedOpBECA");
            ApplyReversedOpBECA(op, register);
        }
        adjoint auto;
        controlled auto;
        controlled adjoint auto;
    }

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Arithmetic.AssertMostSignificantBit".
    operation AssertHighestBit(value : Result, number : LittleEndian) : Unit {
        body (...) {
            _Renamed("Microsoft.Quantum.Arithmetic.AssertHighestBit", "Microsoft.Quantum.Arithmetic.AssertMostSignificantBit");
            AssertMostSignificantBit(value, number);
        }
        adjoint auto;
        controlled auto;
        controlled adjoint auto;
    }

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Arithmetic.LittleEndianAsBigEndian".
    function LittleEndianToBigEndian(input: LittleEndian) : BigEndian {
        _Renamed("LittleEndianToBigEndian", "LittleEndianAsBigEndian");
        return LittleEndianAsBigEndian(input);
    }

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Arithmetic.BigEndianAsLittleEndian".
    function BigEndianToLittleEndian(input: BigEndian) : LittleEndian {
        _Renamed("BigEndianToLittleEndian", "BigEndianAsLittleEndian");
        return BigEndianAsLittleEndian(input);
    }

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Arithmetic.IncrementPhaseByInteger".
    operation IntegerIncrementPhaseLE (increment : Int, target : PhaseLittleEndian) : Unit {
        body (...) {
            _Renamed("Microsoft.Quantum.Arithmetic.IntegerIncrementPhaseLE", "Microsoft.Quantum.Arithmetic.IncrementPhaseByInteger");
            IncrementPhaseByInteger(increment, target);
        }

        adjoint invert;
        controlled distribute;
        controlled adjoint distribute;
    }

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Arithmetic.IncrementByInteger".
    operation IntegerIncrementLE(increment : Int, target : LittleEndian) : Unit {
        body (...) {
            _Renamed("Microsoft.Quantum.Arithmetic.IntegerIncrementLE", "Microsoft.Quantum.Arithmetic.IncrementByInteger");
            IncrementByInteger(increment, target);
        }

        adjoint invert;
        controlled distribute;
        controlled adjoint distribute;
    }

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Measurement.ApplyXorInPlace".
    operation InPlaceXorLE(value : Int, target : LittleEndian) : Unit {
        body (...) {
            _Renamed("Microsoft.Quantum.Canon.InPlaceXorLE", "Microsoft.Quantum.Arithmetic.ApplyXorInPlace");
            ApplyXorInPlace(value, target);
        }
        adjoint auto;
        controlled auto;
        controlled adjoint auto;
    }

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Arithmetic.CompareUsingRippleCarry".
    operation ApplyRippleCarryComparatorLE(x: LittleEndian, y: LittleEndian, output: Qubit) : Unit {
        body (...) {
            _Renamed("Microsoft.Quantum.Canon.ApplyRippleCarryComparatorLE", "Microsoft.Quantum.Arithmetic.CompareUsingRippleCarry");
            CompareUsingRippleCarry(x, y, output);
        }
        adjoint auto;
        controlled auto;
        controlled adjoint auto;
    }

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Arithmetic.ModularIncrementByInteger".
    operation ModularIncrementLE (increment : Int, modulus : Int, target : LittleEndian) : Unit {
        body (...) {
            _Renamed("Microsoft.Quantum.Canon.ModularIncrementLE", "Microsoft.Quantum.Arithmetic.IncrementByModularInteger");
            IncrementByModularInteger(increment, modulus, target);
        }
        adjoint auto;
        controlled auto;
        controlled adjoint auto;
    }

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Arithmetic.IncrementPhaseByModularInteger".
    operation ModularIncrementPhaseLE (increment : Int, modulus : Int, target : PhaseLittleEndian) : Unit {
        body (...) {
            _Renamed("Microsoft.Quantum.Canon.ModularIncrementPhaseLE", "Microsoft.Quantum.Arithmetic.IncrementPhaseByModularInteger");
            IncrementPhaseByModularInteger(increment, modulus, target);
        }
        adjoint auto;
        controlled auto;
        controlled adjoint auto;
    }

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Arithmetic.MultiplyAndAddByModularInteger".
    operation ModularAddProductLE (constMultiplier : Int, modulus : Int, multiplier : LittleEndian, summand : LittleEndian) : Unit {
        body (...) {
            _Renamed("Microsoft.Quantum.Canon.ModularIncrementPhaseLE", "Microsoft.Quantum.Arithmetic.MultiplyAndAddByModularInteger");
            MultiplyAndAddByModularInteger(constMultiplier, modulus, multiplier, summand);
        }
        adjoint auto;
        controlled auto;
        controlled adjoint auto;
    }

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Arithmetic.MultiplyAndAddPhaseByModularInteger".
    operation ModularAddProductPhaseLE (constMultiplier : Int, modulus : Int, multiplier : LittleEndian, phaseSummand : PhaseLittleEndian) : Unit {
        body (...) {
            _Renamed("Microsoft.Quantum.Canon.ModularAddProductPhaseLE", "Microsoft.Quantum.Arithmetic.MultiplyAndAddPhaseByModularInteger");
            MultiplyAndAddPhaseByModularInteger(constMultiplier, modulus, multiplier, phaseSummand);
        }
        adjoint auto;
        controlled auto;
        controlled adjoint auto;
    }

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Arithmetic.MultiplyByModularInteger".
    operation ModularMultiplyByConstantLE(constMultiplier : Int, modulus : Int, multiplier : LittleEndian) : Unit {
        body (...) {
            _Renamed("Microsoft.Quantum.Canon.ModularMultiplyByConstantLE", "Microsoft.Quantum.Arithmetic.MultiplyByModularInteger");
            MultiplyByModularInteger(constMultiplier, modulus, multiplier);
        }
        adjoint auto;
        controlled auto;
        controlled adjoint auto;
    }

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Arithmetic.AssertPhaseLessThan".
    operation AssertLessThanPhaseLE(value : Int, number : PhaseLittleEndian) : Unit {
        body (...) {
            _Renamed("Microsoft.Quantum.Canon.AssertPhaseLessThan", "Microsoft.Quantum.Arithmetic.AssertPhaseLessThan");
            AssertPhaseLessThan(value, number);
        }
        adjoint auto;
        controlled auto;
        adjoint controlled auto;
    }

    // #region Removal of "BE" suffix
    // NB: with Q# 0.6, all arithmetic functionality has been normalized to
    //     take exclusively LittleEndian inputs.

    /// # Deprecated
    /// This operation has been removed.
    operation InPlaceXorBE(value : Int, target : BigEndian) : Unit {
        body (...) {
            _Removed("Microsoft.Quantum.Canon.InPlaceXorBE", "ApplyReversedOpLECA(ApplyXorInPlace(value, _), target)");
            ApplyReversedOpLECA(ApplyXorInPlace(value, _), target);
        }
        adjoint auto;
        controlled auto;
        adjoint controlled auto;
    }

    /// # Deprecated
    /// This operation has been removed.
    operation ApplyRippleCarryComparatorBE(x : BigEndian, y : BigEndian, output : Qubit) : Unit {
        body (...) {
            _Removed(
                "Microsoft.Quantum.Canon.ApplyRippleCarryComparatorBE",
                "CompareUsingRippleCarry(BigEndianAsLittleEndian(x), BigEndianAsLittleEndian(y), output)"
            );
            CompareUsingRippleCarry(BigEndianAsLittleEndian(x), BigEndianAsLittleEndian(y), output);
        }
        adjoint auto;
        controlled auto;
        adjoint controlled auto;
    }

    /// # Deprecated
    /// This operation has been removed.
    operation MeasureIntegerBE(target : BigEndian) : Int {
        _Removed(
            "Microsoft.Quantum.Canon.MeasureIntegerBE",
            "MeasureInteger(BigEndianAsLittleEndian(target))"
        );
        return MeasureInteger(BigEndianAsLittleEndian(target));
    }

    /// # Deprecated
    /// This operation has been removed.
    operation AssertProbIntBE (stateIndex : Int, prob : Double, qubits : BigEndian, tolerance : Double) : Unit {
        _Removed(
            "Microsoft.Quantum.Canon.AssertProbIntBE",
            "AssertProbInt(stateIndex, prob, BigEndianAsLittleEndian(qubits), tolerance)"
        );
        AssertProbInt(stateIndex, prob, BigEndianAsLittleEndian(qubits), tolerance);
    }

    // #endregion

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Arithmetic.CopyMostSignificantBit".
     operation CopyMostSignificantBitLE(from : LittleEndian, target : Qubit) : Unit {
        body (...) {
            _Renamed("Microsoft.Quantum.Canon.CopyMostSignificantBitLE", "Microsoft.Quantum.Arithmetic.CopyMostSignificantBit");
            CopyMostSignificantBit(from, target);
        }

        adjoint invert;
    }

}
