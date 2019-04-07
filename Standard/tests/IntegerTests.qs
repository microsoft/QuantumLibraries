// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Arithmetic {
    open Microsoft.Quantum.Primitive;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Math;
    open Microsoft.Quantum.Diagnostics;
    open Microsoft.Quantum.Measurement;

    operation IntegerAdderTestHelper( IntegerAdder : ( (LittleEndian, LittleEndian, Qubit) => Unit : Controlled), summand1 : Int, summand2 : Int, numberOfQubits : Int ) : Unit {
        body (...) {
            using (register = Qubit[2*numberOfQubits + 1]) {
                mutable actual_carry = 0;
                mutable actual1 = 0;
                mutable actual2 = 0;
                mutable measured_carry = Zero;
                let summand1LE = LittleEndian(register[0 .. numberOfQubits - 1]);
                let summand2LE = LittleEndian(register[numberOfQubits .. 2*numberOfQubits - 1]);
                let carry = register[2*numberOfQubits];
 
                InPlaceXorLE(summand1, summand1LE);
                InPlaceXorLE(summand2, summand2LE);

                IntegerAdder(summand1LE, summand2LE, carry);
 
                let sum = summand1 + summand2;
                let expected = Modulus(sum, 2^numberOfQubits);
                set actual1 = MeasureInteger(summand1LE);
                EqualityFactI(summand1, actual1, $"Expected {summand1}, got {actual1}");
                set actual2 = MeasureInteger(summand2LE);
                EqualityFactI(expected, actual2, $"Expected {expected}, got {actual2}");
                let expected_carry = (sum / 2^numberOfQubits);
                set measured_carry = MResetZ(carry);
                if (measured_carry == One) {set actual_carry = 1;} else {set actual_carry = 0;}
                EqualityFactI(expected_carry, actual_carry, $"Expected {expected_carry}, got {actual_carry}");

                for (numberOfControls in 1..2) { 
                    using (controls = Qubit[numberOfControls]) {
                        InPlaceXorLE(summand1, summand1LE);
                        InPlaceXorLE(summand2, summand2LE);
                        // controls are |0>, no addition is computed
                        (Controlled IntegerAdder) (controls, (summand1LE, summand2LE, carry));
                        set actual1 = MeasureInteger(summand1LE);
                        EqualityFactI(summand1, actual1, $"Expected {summand1}, got {actual1}");
                        set actual2 = MeasureInteger(summand2LE);
                        EqualityFactI(summand2, actual2, $"Expected {expected}, got {actual2}");
                        set measured_carry = MResetZ(carry);
                        if (measured_carry == One) {set actual_carry = 1;} else {set actual_carry = 0;}
                        EqualityFactI(0, actual_carry, $"Expected {0}, got {actual_carry}");
                        InPlaceXorLE(summand1, summand1LE);
                        InPlaceXorLE(summand2, summand2LE);
                        // now controls are set to |1>, addition is computed
                        ApplyToEach(X, controls);
                        (Controlled IntegerAdder) (controls, (summand1LE, summand2LE, carry));
                        set actual1 = MeasureInteger(summand1LE);
                        EqualityFactI(summand1, actual1, $"Expected {summand1}, got {actual1}");
                        set actual2 = MeasureInteger(summand2LE);
                        EqualityFactI(expected, actual2, $"Expected {expected}, got {actual2}");
                        set measured_carry = MResetZ(carry);
                        if (measured_carry == One) {set actual_carry = 1;} else {set actual_carry = 0;}
                        EqualityFactI(expected_carry, actual_carry, $"Expected {expected_carry}, got {actual_carry}");
                        ResetAll(controls);
                    }
                }
            }
        }
    }
 
    operation IntegerAdderExhaustiveTestHelper (IntegerAdder : ( (LittleEndian, LittleEndian, Qubit) => Unit : Controlled), numberOfQubits : Int) : Unit {
        for( summand1 in 0 .. 2^numberOfQubits - 1 ) {
            for( summand2 in 0 .. 2^numberOfQubits - 1 ) {
                IntegerAdderTestHelper(IntegerAdder, summand1, summand2, numberOfQubits);
            }
        }
    }

    operation RippleCarryAdderDTest () : Unit {
        let numberOfQubits = 7;
        let summand1 = 127;
        let summand2 = 17;
        IntegerAdderTestHelper(RippleCarryAdderD, summand1, summand2, numberOfQubits);
    }

    operation RippleCarryAdderDExhaustiveTestReversible () : Unit {
        for (numberOfQubits in 3..6) {
            IntegerAdderExhaustiveTestHelper (RippleCarryAdderD, numberOfQubits);
        }
    }

    operation RippleCarryAdderDTestReversible () : Unit {
        let numberOfQubits = 20;
        let summand1 = 823709; 
        let summand2 = 88487;  
        IntegerAdderTestHelper(RippleCarryAdderD, summand1, summand2, numberOfQubits);
    }

    operation RippleCarryAdderCDKMExhaustiveTest () : Unit {
        let numberOfQubits = 4;
        IntegerAdderExhaustiveTestHelper (RippleCarryAdderCDKM, numberOfQubits);    
    }

    operation RippleCarryAdderCDKMTestReversible () : Unit {
        let numberOfQubits = 20;
        let summand1 = 823709; 
        let summand2 = 88487;  
        IntegerAdderTestHelper(RippleCarryAdderCDKM, summand1, summand2, numberOfQubits);
    }

    operation RippleCarryAdderCDKMExhaustiveTestReversible () : Unit {
        for (numberOfQubits in 3..6) {
            IntegerAdderExhaustiveTestHelper (RippleCarryAdderCDKM, numberOfQubits); 
        }
    }

    operation RippleCarryAdderTTKExhaustiveTest () : Unit {
        let numberOfQubits = 4;
        IntegerAdderExhaustiveTestHelper (RippleCarryAdderTTK, numberOfQubits);
    }

    operation RippleCarryAdderTTKExhaustiveTestReversible () : Unit {
        for (numberOfQubits in 1..6){
            IntegerAdderExhaustiveTestHelper (RippleCarryAdderTTK, numberOfQubits);
        }
    }


    operation IntegerAdderNoCarryTestHelper( IntegerAdder : ( (LittleEndian, LittleEndian) => Unit : Controlled), summand1 : Int, summand2 : Int, numberOfQubits : Int ) : Unit {
        using (register = Qubit[2*numberOfQubits]) {
            mutable actual1 = 0;
            mutable actual2 = 0;
            let summand1LE = LittleEndian(register[0 .. numberOfQubits - 1]);
            let summand2LE = LittleEndian(register[numberOfQubits .. 2*numberOfQubits - 1]);
 
            InPlaceXorLE(summand1, summand1LE);
            InPlaceXorLE(summand2, summand2LE);

            IntegerAdder(summand1LE, summand2LE);
 
            let sum = summand1 + summand2;
            let expected = Modulus(sum, 2^numberOfQubits);
            set actual1 = MeasureInteger(summand1LE);
            EqualityFactI(summand1, actual1, $"Expected {summand1}, got {actual1}");
            set actual2 = MeasureInteger(summand2LE);
            EqualityFactI(expected, actual2, $"Expected {expected}, got {actual2}");
            let expected_carry = (sum / 2^numberOfQubits);

            for (numberOfControls in 1..2) { 
                using (controls = Qubit[numberOfControls]) {
                    InPlaceXorLE(summand1, summand1LE);
                    InPlaceXorLE(summand2, summand2LE);
                    // controls are |0>, no addition is computed
                    (Controlled IntegerAdder) (controls, (summand1LE, summand2LE));
                    set actual1 = MeasureInteger(summand1LE);
                    EqualityFactI(summand1, actual1, $"Expected {summand1}, got {actual1}");
                    set actual2 = MeasureInteger(summand2LE);
                    EqualityFactI(summand2, actual2, $"Expected {expected}, got {actual2}");
                    InPlaceXorLE(summand1, summand1LE);
                    InPlaceXorLE(summand2, summand2LE);
                    // now controls are set to |1>, addition is computed
                    ApplyToEach(X, controls);
                    (Controlled IntegerAdder) (controls, (summand1LE, summand2LE));
                    set actual1 = MeasureInteger(summand1LE);
                    EqualityFactI(summand1, actual1, $"Expected {summand1}, got {actual1}");
                    set actual2 = MeasureInteger(summand2LE);
                    EqualityFactI(expected, actual2, $"Expected {expected}, got {actual2}");
                    ResetAll(controls);
                }
            }
        }
    }

    operation IntegerAdderNoCarryExhaustiveTestHelper (IntegerAdder : ( (LittleEndian, LittleEndian) => Unit : Controlled), numberOfQubits : Int) : Unit {
        for( summand1 in 0 .. 2^numberOfQubits - 1 ) {
            for( summand2 in 0 .. 2^numberOfQubits - 1 ) {
                IntegerAdderNoCarryTestHelper(IntegerAdder, summand1, summand2, numberOfQubits);
            }
        }
    }

    operation RippleCarryAdderNoCarryTTKTestReversible () : Unit {
        let numberOfQubits = 10;
        let summand1 = 1021; 
        let summand2 = 973; 
        IntegerAdderNoCarryTestHelper(RippleCarryAdderNoCarryTTK, summand1, summand2, numberOfQubits);
    }

    operation RippleCarryAdderNoCarryTTKExhaustiveTest () : Unit {
        let numberOfQubits = 4;
        IntegerAdderNoCarryExhaustiveTestHelper (RippleCarryAdderNoCarryTTK, numberOfQubits);
    }

    operation RippleCarryAdderNoCarryTTKExhaustiveTestReversible () : Unit {
        for (numberOfQubits in 1..6) {
            IntegerAdderNoCarryExhaustiveTestHelper (RippleCarryAdderNoCarryTTK, numberOfQubits);
        }
    }


    operation GreaterThanTestHelper( integer1 : Int, integer2 : Int, numberOfQubits : Int ) : Unit {
        using (register = Qubit[2*numberOfQubits+1]) {
            mutable actual1 = 0;
            mutable actual2 = 0;
            mutable actualr = Zero;
            mutable gt = Zero;
            let integer1LE = LittleEndian(register[0 .. numberOfQubits - 1]);
            let integer2LE = LittleEndian(register[numberOfQubits .. 2*numberOfQubits - 1]);
            let result = register[2*numberOfQubits];
 
            InPlaceXorLE(integer1, integer1LE);
            InPlaceXorLE(integer2, integer2LE);

            GreaterThan(integer1LE, integer2LE, result);

            if (integer1 > integer2) {set gt = One;} 
            set actual1 = MeasureInteger(integer1LE);
            EqualityFactI(integer1, actual1, $"Expected {integer1}, got {actual1}");
            set actual2 = MeasureInteger(integer2LE);
            EqualityFactI(integer2, actual2, $"Expected {integer2}, got {actual2}");
            set actualr = M(result);
            EqualityFactB((gt == actualr), true, $"Expected {gt}, got {actualr}");

            Reset(result);
            for (numberOfControls in 1..2) { 
                using (controls = Qubit[numberOfControls]) {
                    InPlaceXorLE(integer1, integer1LE);
                    InPlaceXorLE(integer2, integer2LE);
                    (Controlled GreaterThan) (controls, (integer1LE, integer2LE, result));

                    set actual1 = MeasureInteger(integer1LE);
                    EqualityFactI(integer1, actual1, $"Expected {integer1}, got {actual1}");
                    set actual2 = MeasureInteger(integer2LE);
                    EqualityFactI(integer2, actual2, $"Expected {integer2}, got {actual2}");
                    set actualr = M(result);
                    EqualityFactB((actualr == Zero), true, $"Expected Zero, got {actualr}");

                    ApplyToEach(X, controls);
                    InPlaceXorLE(integer1, integer1LE);
                    InPlaceXorLE(integer2, integer2LE);
                    (Controlled GreaterThan) (controls, (integer1LE, integer2LE, result));

                    set actual1 = MeasureInteger(integer1LE);
                    EqualityFactI(integer1, actual1, $"Expected {integer1}, got {actual1}");
                    set actual2 = MeasureInteger(integer2LE);
                    EqualityFactI(integer2, actual2, $"Expected {integer2}, got {actual2}");
                    set actualr = M(result);
                    EqualityFactB((gt == actualr), true, $"Expected {gt}, got {actualr}");

                    ResetAll(controls);
                    Reset(result);
                }
            }
        }
    }

    operation GreaterThanExhaustiveTestReversible () : Unit {
        for (numberOfQubits in 1..5) {
            for (integer1 in 0..2^numberOfQubits-1) {
                for (integer2 in 0..2^numberOfQubits-1) {
                    GreaterThanTestHelper(integer1, integer2, numberOfQubits);
                }
            }
        }
    }

}
