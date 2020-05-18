// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {
    open Microsoft.Quantum.Arithmetic;
    open Microsoft.Quantum.Measurement;

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Arithmetic.ApplyReversedOpLEA".
    @Deprecated("Microsoft.Quantum.Arithmetic.ApplyReversedOpLEA")
    operation ApplyReversedOpLittleEndianA(op : (LittleEndian => Unit is Adj), register : BigEndian) : Unit is Adj {
        ApplyReversedOpLEA(op, register);
    }

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Arithmetic.ApplyReversedOpLEC".
    @Deprecated("Microsoft.Quantum.Arithmetic.ApplyReversedOpLEC")
    operation ApplyReversedOpLittleEndianC(op : (LittleEndian => Unit is Ctl), register : BigEndian) : Unit is Ctl {
        ApplyReversedOpLEC(op, register);
    }

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Arithmetic.ApplyReversedOpLEC".
    @Deprecated("Microsoft.Quantum.Arithmetic.ApplyReversedOpLECA")
    operation ApplyReversedOpLittleEndianCA(op : (LittleEndian => Unit is Adj + Ctl), register : BigEndian) : Unit is Adj + Ctl {
        ApplyReversedOpLECA(op, register);
    }

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Arithmetic.ApplyReversedOpBEA".
    @Deprecated("Microsoft.Quantum.Arithmetic.ApplyReversedOpBEA")
    operation ApplyReversedOpBigEndianA(op : (BigEndian => Unit is Adj), register : LittleEndian) : Unit is Adj {
        ApplyReversedOpBEA(op, register);
    }

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Arithmetic.ApplyReversedOpBEC".
    @Deprecated("Microsoft.Quantum.Arithmetic.ApplyReversedOpBEC")
    operation ApplyReversedOpBigEndianC(op : (BigEndian => Unit is Ctl), register : LittleEndian) : Unit is Ctl {
        ApplyReversedOpBEC(op, register);
    }

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Arithmetic.ApplyReversedOpBEC".
    @Deprecated("Microsoft.Quantum.Arithmetic.ApplyReversedOpBECA")
    operation ApplyReversedOpBigEndianCA(op : (BigEndian => Unit is Adj + Ctl), register : LittleEndian) : Unit is Ctl + Adj {
        ApplyReversedOpBECA(op, register);
    }

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Arithmetic.AssertMostSignificantBit".
    @Deprecated("Microsoft.Quantum.Arithmetic.AssertMostSignificantBit")
    operation AssertHighestBit(value : Result, number : LittleEndian) : Unit is Ctl + Adj {
        AssertMostSignificantBit(value, number);
    }

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Arithmetic.LittleEndianAsBigEndian".
    @Deprecated("Microsoft.Quantum.Arithmetic.LittleEndianAsBigEndian")
    function LittleEndianToBigEndian(input: LittleEndian) : BigEndian {
        return LittleEndianAsBigEndian(input);
    }

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Arithmetic.BigEndianAsLittleEndian".
    @Deprecated("Microsoft.Quantum.Arithmetic.BigEndianAsLittleEndian")
    function BigEndianToLittleEndian(input: BigEndian) : LittleEndian {
        return BigEndianAsLittleEndian(input);
    }

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Arithmetic.IncrementPhaseByInteger".
    @Deprecated("Microsoft.Quantum.Arithmetic.IncrementPhaseByInteger")
    operation IntegerIncrementPhaseLE (increment : Int, target : PhaseLittleEndian) : Unit is Adj + Ctl {
        IncrementPhaseByInteger(increment, target);
    }

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Arithmetic.IncrementByInteger".
    @Deprecated("Microsoft.Quantum.Arithmetic.IncrementByInteger")
    operation IntegerIncrementLE(increment : Int, target : LittleEndian) : Unit is Adj + Ctl {
        IncrementByInteger(increment, target);
    }

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Measurement.ApplyXorInPlace".
    @Deprecated("Microsoft.Quantum.Measurement.ApplyXorInPlace")
    operation InPlaceXorLE(value : Int, target : LittleEndian) : Unit is Adj + Ctl {
        ApplyXorInPlace(value, target);
    }

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Arithmetic.CompareUsingRippleCarry".
    @Deprecated("Microsoft.Quantum.Arithmetic.CompareUsingRippleCarry")
    operation ApplyRippleCarryComparatorLE(x: LittleEndian, y: LittleEndian, output: Qubit) : Unit is Adj + Ctl {
        CompareUsingRippleCarry(x, y, output);
    }

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Arithmetic.ModularIncrementByInteger".
    @Deprecated("Microsoft.Quantum.Arithmetic.ModularIncrementByInteger")
    operation ModularIncrementLE (increment : Int, modulus : Int, target : LittleEndian) : Unit is Adj + Ctl {
        IncrementByModularInteger(increment, modulus, target);
    }

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Arithmetic.IncrementPhaseByModularInteger".
    @Deprecated("Microsoft.Quantum.Arithmetic.IncrementPhaseByModularInteger")
    operation ModularIncrementPhaseLE (increment : Int, modulus : Int, target : PhaseLittleEndian) : Unit is Adj + Ctl {
        IncrementPhaseByModularInteger(increment, modulus, target);
    }

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Arithmetic.MultiplyAndAddByModularInteger".
    @Deprecated("Microsoft.Quantum.Arithmetic.MultiplyAndAddByModularInteger")
    operation ModularAddProductLE (constMultiplier : Int, modulus : Int, multiplier : LittleEndian, summand : LittleEndian) : Unit is Adj + Ctl {
        MultiplyAndAddByModularInteger(constMultiplier, modulus, multiplier, summand);
    }

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Arithmetic.MultiplyAndAddPhaseByModularInteger".
    @Deprecated("Microsoft.Quantum.Arithmetic.MultiplyAndAddPhaseByModularInteger")
    operation ModularAddProductPhaseLE (constMultiplier : Int, modulus : Int, multiplier : LittleEndian, phaseSummand : PhaseLittleEndian) : Unit is Adj + Ctl {
        MultiplyAndAddPhaseByModularInteger(constMultiplier, modulus, multiplier, phaseSummand);
    }

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Arithmetic.MultiplyByModularInteger".
    @Deprecated("Microsoft.Quantum.Arithmetic.MultiplyByModularInteger")
    operation ModularMultiplyByConstantLE(constMultiplier : Int, modulus : Int, multiplier : LittleEndian) : Unit is Adj + Ctl {
        MultiplyByModularInteger(constMultiplier, modulus, multiplier);
    }

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Arithmetic.AssertPhaseLessThan".
    @Deprecated("Microsoft.Quantum.Arithmetic.AssertPhaseLessThan")
    operation AssertLessThanPhaseLE(value : Int, number : PhaseLittleEndian) : Unit is Adj + Ctl {
        AssertPhaseLessThan(value, number);
    }

    /// # Deprecated
    /// This operation has been removed.
    @Deprecated("ApplyReversedOpLECA(ApplyXorInPlace(value, _), target)")
    operation InPlaceXorBE(value : Int, target : BigEndian) : Unit is Adj + Ctl {
        ApplyReversedOpLECA(ApplyXorInPlace(value, _), target);
    }

    /// # Deprecated
    /// Please use @"microsoft.quantum.canon.applycnotchain".
    @Deprecated("Microsoft.Quantum.Canon.ApplyCNOTChain")
    operation CascadeCNOT (register : Qubit[]) : Unit is Adj + Ctl {
        Microsoft.Quantum.Canon.ApplyCNOTChain(register);
    }

    /// # Deprecated
    /// Please use @"microsoft.quantum.arithmetic.applymajorityinplace".
    @Deprecated("Microsoft.Quantum.Arithmetic.ApplyMajorityInPlace")
    operation InPlaceMajority(output: Qubit, input: Qubit[])
    : Unit is Adj + Ctl {
        ApplyMajorityInPlace(output, input);
    }

    @Deprecated("Microsoft.Quantum.Canon.ApplyCCNOTChain")
    operation CascadeCCNOT(register : Qubit[], targets : Qubit[])
    : Unit is Adj + Ctl {
        ApplyCCNOTChain(register, targets);
    }

}
