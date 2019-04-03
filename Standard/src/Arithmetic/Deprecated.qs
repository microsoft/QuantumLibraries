// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {
    open Microsoft.Quantum.Arithmetic;
    open Microsoft.Quantum.Measurement;

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Arithmetic.ApplyReversedOpLE".
    operation ApplyReversedOpLittleEndian(op : (LittleEndian => Unit), register : BigEndian) : Unit {
        Renamed("ApplyReversedOpLittleEndian", "ApplyReversedOpLE");
        ApplyReversedOpLE(op, register);
    }

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Arithmetic.ApplyReversedOpLEA".
    operation ApplyReversedOpLittleEndianA(op : (LittleEndian => Unit : Adjoint), register : BigEndian) : Unit {
        body (...) {
            Renamed("ApplyReversedOpLittleEndianA", "ApplyReversedOpLEA");
            ApplyReversedOpLEA(op, register);
        }
        adjoint auto;
    }

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Arithmetic.ApplyReversedOpLEC".
    operation ApplyReversedOpLittleEndianC(op : (LittleEndian => Unit : Controlled), register : BigEndian) : Unit {
        body (...) {
            Renamed("ApplyReversedOpLittleEndianC", "ApplyReversedOpLEC");
            ApplyReversedOpLEC(op, register);
        }
        controlled auto;
    }

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Arithmetic.ApplyReversedOpLEC".
    operation ApplyReversedOpLittleEndianCA(op : (LittleEndian => Unit : Adjoint, Controlled), register : BigEndian) : Unit {
        body (...) {
            Renamed("ApplyReversedOpLittleEndianCA", "ApplyReversedOpLECA");
            ApplyReversedOpLECA(op, register);
        }
        adjoint auto;
        controlled auto;
        controlled adjoint auto;
    }

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Arithmetic.ApplyReversedOpBE".
    operation ApplyReversedOpBigEndian(op : (BigEndian => Unit), register : LittleEndian) : Unit {
        Renamed("ApplyReversedOpBigEndian", "ApplyReversedOpBE");
        ApplyReversedOpBE(op, register);
    }

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Arithmetic.ApplyReversedOpBEA".
    operation ApplyReversedOpBigEndianA(op : (BigEndian => Unit : Adjoint), register : LittleEndian) : Unit {
        body (...) {
            Renamed("ApplyReversedOpBigEndianA", "ApplyReversedOpBEA");
            ApplyReversedOpBEA(op, register);
        }
        adjoint auto;
    }
    
    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Arithmetic.ApplyReversedOpBEC".
    operation ApplyReversedOpBigEndianC(op : (BigEndian => Unit : Controlled), register : LittleEndian) : Unit {
        body (...) {
            Renamed("ApplyReversedOpBigEndianC", "ApplyReversedOpBEC");
            ApplyReversedOpBEC(op, register);
        }
        controlled auto;
    }

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Arithmetic.ApplyReversedOpBEC".
    operation ApplyReversedOpBigEndianCA(op : (BigEndian => Unit : Adjoint, Controlled), register : LittleEndian) : Unit {
        body (...) {
            Renamed("ApplyReversedOpBigEndianCA", "ApplyReversedOpBECA");
            ApplyReversedOpBECA(op, register);
        }
        adjoint auto;
        controlled auto;
        controlled adjoint auto;
    }

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Arithmetic.AssertProbIntLE".
    operation AssertProbInt(stateIndex : Int, expected : Double, qubits : LittleEndian, tolerance : Double) : Unit {
        Renamed("Microsoft.Quantum.Arithmetic.AssertProbInt", "Microsoft.Quantum.Arithmetic.AssertProbIntLE");
        AssertProbIntLE(stateIndex, expected, qubits, tolerance);
    }

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Arithmetic.AssertMostSignificantBitLE".
    operation AssertHighestBit(value : Result, number : LittleEndian) : Unit {
        body (...) {
            Renamed("Microsoft.Quantum.Arithmetic.AssertHighestBit", "Microsoft.Quantum.Arithmetic.AssertMostSignificantBitLE");
            AssertMostSignificantBitLE(value, number);
        }
        adjoint auto;
        controlled auto;
        controlled adjoint auto;
    }

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Arithmetic.LittleEndianAsBigEndian".
    function LittleEndianToBigEndian(input: LittleEndian) : BigEndian {
        Renamed("LittleEndianToBigEndian", "LittleEndianAsBigEndian");
        return LittleEndianAsBigEndian(input);
    }

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Arithmetic.BigEndianAsLittleEndian".
    function BigEndianToLittleEndian(input: BigEndian) : LittleEndian {
        Renamed("BigEndianToLittleEndian", "BigEndianAsLittleEndian");
        return BigEndianAsLittleEndian(input);
    }

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Arithmetic.IncrementPhaseByInteger".
    operation IntegerIncrementPhaseLE (increment : Int, target : PhaseLittleEndian) : Unit {
        body (...) {
            Renamed("Microsoft.Quantum.Arithmetic.IntegerIncrementPhaseLE", "Microsoft.Quantum.Arithmetic.IncrementPhaseByInteger");
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
            Renamed("Microsoft.Quantum.Arithmetic.IntegerIncrementLE", "Microsoft.Quantum.Arithmetic.IncrementByInteger");
            IncrementByInteger(increment, target);
        }

        adjoint invert;
        controlled distribute;
        controlled adjoint distribute;
    }

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Measurement.MeasureIntegerLE".
    operation MeasureInteger(target : LittleEndian) : Int {
        Renamed("Microsoft.Quantum.Canon.MeasureInteger", "Microsoft.Quantum.Measurement.MeasureIntegerLE");
        return MeasureIntegerLE(target);
    }

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Measurement.ApplyXorInPlace".
    operation InPlaceXorLE(value : Int, target : LittleEndian) : Unit {
        body (...) {
            Renamed("Microsoft.Quantum.Canon.InPlaceXorLE", "Microsoft.Quantum.Arithmetic.ApplyXorInPlace");
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
            Renamed("Microsoft.Quantum.Canon.ApplyRippleCarryComparatorLE", "Microsoft.Quantum.Arithmetic.CompareUsingRippleCarry");
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
            Renamed("Microsoft.Quantum.Canon.ModularIncrementLE", "Microsoft.Quantum.Arithmetic.IncrementByModularInteger");
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
            Renamed("Microsoft.Quantum.Canon.ModularIncrementPhaseLE", "Microsoft.Quantum.Arithmetic.IncrementPhaseByModularInteger");
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
            Renamed("Microsoft.Quantum.Canon.ModularIncrementPhaseLE", "Microsoft.Quantum.Arithmetic.MultiplyAndAddByModularInteger");
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
            Renamed("Microsoft.Quantum.Canon.ModularAddProductPhaseLE", "Microsoft.Quantum.Arithmetic.MultiplyAndAddPhaseByModularInteger");
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
            Renamed("Microsoft.Quantum.Canon.ModularMultiplyByConstantLE", "Microsoft.Quantum.Arithmetic.MultiplyByModularInteger");
            MultiplyByModularInteger(constMultiplier, modulus, multiplier);
        }
        adjoint auto;
        controlled auto;
        controlled adjoint auto;
    }


    // #region Removal of "BE" suffix
    // NB: with Q# 0.6, all arithmetic functionality has been normalized to
    //     take exclusively LittleEndian inputs.

    /// # Deprecated
    /// This operation has been removed.
    operation InPlaceXorBE(value : Int, target : BigEndian) : Unit {
        body (...) {
            Removed("Microsoft.Quantum.Canon.InPlaceXorBE", "ApplyReversedOpLECA(ApplyXorInPlace(value, _), target)");
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
            Removed(
                "Microsoft.Quantum.Canon.ApplyRippleCarryComparatorBE",
                "ApplyRippleCarryComparatorLE(BigEndianAsLittleEndian(x), BigEndianAsLittleEndian(y), output)"
            );
            ApplyRippleCarryComparatorLE(BigEndianAsLittleEndian(x), BigEndianAsLittleEndian(y), output);
        }
        adjoint auto;
        controlled auto;
        adjoint controlled auto;
    }

    // #endregion

}
