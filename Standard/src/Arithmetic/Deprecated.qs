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
    /// Please use @"Microsoft.Quantum.Arithmetic.IncrementByIntegerPhaseLE".
    operation IntegerIncrementPhaseLE (increment : Int, target : PhaseLittleEndian) : Unit {
        body (...) {
            Renamed("Microsoft.Quantum.Arithmetic.IntegerIncrementPhaseLE", "Microsoft.Quantum.Arithmetic.IncrementByIntegerPhaseLE");
            IncrementByIntegerPhaseLE(increment, target);
        }

        adjoint invert;
        controlled distribute;
        controlled adjoint distribute;
    }

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Arithmetic.IncrementByIntegerLE".
    operation IntegerIncrementLE(increment : Int, target : LittleEndian) : Unit {
        body (...) {
            Renamed("Microsoft.Quantum.Arithmetic.IntegerIncrementLE", "Microsoft.Quantum.Arithmetic.IncrementByIntegerLE");
            IncrementByIntegerLE(increment, target);
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

}
