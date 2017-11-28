namespace Microsoft.Quantum.Tests {
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Primitive;

    operation PrepareIntegerTest() : () {
        body {
            using (register = Qubit[4]) {
                PrepareInteger(4, register);
                AssertIntEqual(4, MeasureInteger(register), "Did not measure the integer we expected.");
                ApplyToEach(Reset, register);
            }
        }
    }
    
    operation MeasureIntegerTest() : () {
        body {
            using (register = Qubit[4]) { 
                PrepareInteger(4, register); 
                let result = MeasureInteger(register); 
                AssertIntEqual(result, 4, "Did not measure the integer we expected.");
                ApplyToEach(Reset, register);
            )
        }
    }
}
