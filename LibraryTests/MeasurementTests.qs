namespace Microsoft.Quantum.Tests {
    open Microsoft.Quantum.Primitive;
	open Microsoft.Quantum.Samples.Measurement;

    operation MeasurementTest () : () {
        body { 
            mutable resOne = Zero;
            set resOne = MeasurementOneQubit();
			mutable resTwo = (Zero, Zero);
            set resTwo = MeasurementTwoQubits();
			mutable resBell = (Zero, Zero);
            set resBell = MeasurementBellBasis();			
        }
    }
}