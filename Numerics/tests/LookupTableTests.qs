namespace Microsoft.Quantum.Tests {
    open Microsoft.Quantum.Arithmetic;
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Math;
    open Microsoft.Quantum.Convert;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Diagnostics;
    open Microsoft.Quantum.Random;

    /// # Summary
    /// Example of exponential operation with [0, 9.23] with input-eps 1e-3, output-eps 1e-3
    /// Input should require something like 5 integer bits and 10 fractional bits
    /// Output should require something like 15 integer bits and 10 fractional bits
    operation ExponentialExample(numSwapBits: Int, input : FixedPoint, output : FixedPoint) : Unit {
        let func = ExpD;
        let xMin = 0.0;
        let xMax = 9.23;
        let epsIn = 1e-3;
        let epsOut = 1e-3;

        let lookup = ApplyFunctionWithLookupTable(ExpD, (xMin, xMax), epsIn, epsOut, numSwapBits);
        // Check that input and output registers are the expected size
        EqualityFactI(lookup::IntegerBitsIn + lookup::FractionalBitsIn, 15, "Number of input bits is incorrect");
        EqualityFactI(lookup::IntegerBitsOut + lookup::FractionalBitsOut, 25, "Number of output bits is incorrect");
        lookup::Apply(input, output);
    }
    
    /// # Summary
    /// Example to call Exponential with pre-computed number of input and
    /// output qubits from ExponentialExample
    operation EstimateExponentialInstance(numSwapBits : Int) : Unit {
        use input = Qubit[15];
        use output = Qubit[25];
        let inputFxP = FixedPoint(5, input);
        let outputFxP = FixedPoint(15, output);

        ExponentialExample(numSwapBits, inputFxP, outputFxP);
    }

    //Tests for some examples
    @Test("ToffoliSimulator")
    operation TestExponentialExample() : Unit {
        use input = Qubit[15];
        use output = Qubit[25];
        let inputFxP = FixedPoint(5, input);
        let outputFxP = FixedPoint(15, output);

        for i in 0..9 {
            let inputValue = DrawRandomDouble(9.23*IntAsDouble(i)/10.0, 9.23*IntAsDouble(i+1)/10.0);

            PrepareFxP(inputValue, inputFxP);
            ExponentialExample(5, inputFxP, outputFxP);
            
            let inResult = MeasureFxP(inputFxP);
            let outResult = MeasureFxP(outputFxP);

            let expected = ExpD(inResult);

            EqualityWithinToleranceFact(inResult, inputValue, 1e-3);
            EqualityWithinToleranceFact(outResult, expected, 1e-3);

        }
    }

    @Test("ToffoliSimulator")
    operation TestFxPRegisterValues(): Unit {
        let func = ExpD;
        let xMin = -1.0;
        let xMax = 8.125;
        let epsIn = 0.125;
        let epsOut = 0.25;
 
        let lookup = ApplyFunctionWithLookupTable(func, (xMin, xMax), epsIn, epsOut, 0);

        use inputRegister = Qubit[lookup::IntegerBitsIn+lookup::FractionalBitsIn];
        use outputRegister = Qubit[lookup::IntegerBitsOut+lookup::FractionalBitsOut];

        let inputFxP = FixedPoint(lookup::IntegerBitsIn, inputRegister);
        let outputFxP = FixedPoint(lookup::IntegerBitsOut, outputRegister);

        let samplePoints = 20;

        for i in 0..samplePoints {
            let x = xMin + IntAsDouble(i) * (xMax - xMin) / IntAsDouble(samplePoints);
            PrepareFxP(x, inputFxP);
            lookup::Apply(inputFxP, outputFxP);
            let inputResult = MeasureFxP(inputFxP);
            let outputResult = MeasureFxP(outputFxP);

            EqualityWithinToleranceFact(inputResult, x, epsIn);
            EqualityWithinToleranceFact(outputResult, func(inputResult), epsOut);

        }
    }
}
