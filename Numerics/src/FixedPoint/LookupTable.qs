namespace Microsoft.Quantum.Arithmetic {
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Convert;
    open Microsoft.Quantum.Diagnostics;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Math;

    /// # Summary
    /// The return type when making a lookup table. This contains the operation that 
    /// makes the lookup table circuit, as well as all the parameters required to make 
    /// the two FixedPoint registers that need to be used as inputs and outputs to the 
    /// operator. 
    ///
    /// # Remarks
    /// The reason we have this return type structure is so that the operator is similar 
    /// to the other typical Q# arithmetic function implementations (a larger discussion
    /// can had as to whether that can be changed)
    newtype FunctionWithLookupTable = ( 
        IntegerBitsIn: Int,
        FractionalBitsIn: Int,
        IntegerBitsOut: Int,
        FractionalBitsOut: Int,
        Apply: (FixedPoint, FixedPoint) => Unit is Adj
    );

    /// # Summary
    /// This function creates a select-swap lookup table operator for the function that you want to approximate, as well as
    /// the parameters required to make the two FixedPoint registers that need to be used as inputs to the operator. 
    /// This is so that it is in similar to the other typical Q# arithmetic function (a larger discussion can be had 
    /// as to whether that can be changed). The circuit for the operator can be found in Fig. 1c in arXiv:1812.00954.
    /// 
    /// # Remarks
    /// The operator guarantees that given an input value $x$ and a function $f(x)$, 
    /// it will compute $\hat{f}(\hat{x})$ where $\hat{f}$ is an approximation of $f$ with a maximum error of epsOut and $\hat{x}$ is an
    /// approximation of the input value $\hat{x}$ with a maximum error of epsIn. This is useful for most reasonably behaved 
    /// functions, but note that it computes $\hat{f}(\hat{x})$ and not $\hat{f}(x)$ so if the domain function is very oscillatory and/or
    /// has funky derivatives then it may have high errors.
    ///
    /// # Input
    /// ## func
    /// The Q# arithmetic function that you want to implement with the lookup table
    /// ## domain
    /// A tuple consisting of the minimum and maximum values of the input values to the function
    /// ## epsIn
    /// The maximum allowed error of the input value to the computation (i.e. |x'-x|)
    /// ## epsOut
    /// The maximum allowed error of the output without taking into account the error in input value (i.e. |f'(x')-f(x')|)
    /// ## numSwapBits
    /// The number of bits of the input register that will be used in the SWAP section of the circuits. Another way of looking
    /// at this is that in step in the SELECT section of the circuit in Fig 1c of arXiv:1812.00954, we will encode 2^numSwapBits
    /// encoded 
    ///
    /// # Example
    /// The following code creates a quantum operation based on `ExpD` in the (inclusive) range from `-5.0` to `5.0` with an input error of `1e-3` and an output error of `1e-4`. It uses `2` SWAP bits for the implementation.
    ///
    /// ```qsharp
    /// // Create operation from lookup table
    /// let domain = (-5.0, 5.0);
    /// let epsIn = 1e-3;
    /// let epsOut = 1e-4;
    ///
    /// let lookup = ApplyFunctionWithLookupTable(ExpD, domain, epsIn, epsOut, 2);
    ///
    /// // Allocate qubits
    /// use input = Qubit[lookup::IntegerBitsIn + lookup::FractionalBitsIn];
    /// use output = Qubit[lookup::IntegerBitsOut + lookup::FractionalBitsOut];
    ///
    /// // Represent qubit registers as fixed points
    /// let inputFxP = FixedPoint(lookup::IntegerBitsIn, input);
    /// let outputFxP = FixedPoint(lookup::IntegerBitsOut, output);
    ///
    /// // Apply operation
    /// lookup::Apply(inputFxP, outputFxP);
    /// ```
    function ApplyFunctionWithLookupTable(func: Double -> Double, domain: (Double, Double), epsIn: Double, epsOut: Double, numSwapBits: Int): FunctionWithLookupTable {

        // First step is to find the number of integer bits (pIn) and fractional bits (qIn) required for the input based on the
        // domain and error tolerance (espIn). To find the value of pIn, we have to check both the
        // lower and upper bound of the domain to see which one requires more bits, then assign the larger one as pIn.
        // To find qIn we compute minimum number of fractional bits required to represent epsIn.
        let (minIn, maxIn) = domain;

        let pLower = BitSizeI(Ceiling(AbsD(minIn)));
        let pUpper = BitSizeI(Ceiling(AbsD(maxIn)));
        let pIn = MaxI(pLower, pUpper) + 1; // The +1 is for the sign bit

        let qIn = Ceiling(Lg(1.0/epsIn));

        // We have now computed the number of integer and fractional bits required for the input of our lookup table. Next we compute
        // The output number of integer and fractional bits required. For the number of fractional bits (qOut), we can
        // simply look at the minimum number of bits required to represent epsOut
        let qOut = Ceiling(Lg(1.0/epsOut));


        // For the number of integer bits required for the output, we have to iterate through all the possible values of the function
        // and find the one with the largest absolute value. For that we first create the fixed point approximations of minIn and maxIn
        // given epsIn (using the previously computed pIn and qIn). Then we compute how many different input values (numValues) are there between
        // minIn and maxIn (given our number of input qubits). And finally we evaluate the function at all those values to get the number with
        // the largest absolute value

        // Compute approximations of minIn and maxIn
        let minInFxP = DoubleAsFixedPoint(pIn, qIn, minIn);
        let maxInFxP = DoubleAsFixedPoint(pIn, qIn, maxIn);

        // Compute number of values in between minIn and maxIn
        let deltaIn = 1.0/PowD(2.0, IntAsDouble(qIn));
        let numValues = Truncate((maxInFxP - minInFxP) / deltaIn) + 1;

        // Go through each value, compute the number of integer bits required, and update pOut if it's bigger than
        // current pOut. We also store the output values since we will be using them when creating the output part of the
        // lookup table circuit
        mutable outValues = [0.0, size=numValues]; // List to store all the output values (initialized at all 0s)
        mutable inValueFxP = minInFxP; // Starting input value
        mutable pOut = 0; // Set initial pOut value which will be updated in loop below
        for i in 0..numValues-1 {
            // First a quick check to see that the enumaration is going correctly, i.e. that we are hitting all the values in order
            let inAddress = BoolArrayAsInt(FixedPointAsBoolArray(pIn, qIn, inValueFxP - minInFxP));
            EqualityFactI(inAddress, i, $"Unexpected address in enumeration");

            // Now we compute the output value, compute the number of integer bits it has and see if it is bigger than our current pOut
            let outValue = func(inValueFxP);
            set outValues w/= i <- outValue; // this is the syntax to say "outValues = outValues but with the ith index as outValue"
            set pOut = MaxI(pOut, BitSizeI(Ceiling(AbsD(outValue)))+1); //the +1 is for the sign bit
            set inValueFxP += deltaIn;
        }

        //So we have now computed the number of integer bits for the output values. Now all that's left is to make the circuit!

        // We first create a list of FixedPoints with all the outValues
        let outValuesFxP = Mapped(DoubleAsFixedPoint(pOut, qOut, _), outValues);

        // Next we map outValuesFP to bitstrings
        let outBits = Mapped(FixedPointAsBoolArray(pOut, qOut, _), outValues);
        // Message($"{outBits}");

        // Now we use the fixed point approximation of the minimum value of the input
        // and the list of output bit values to make the operation lookupOperation: (FixedPoint, FixedPoint) => Unit
        // More comments on how that's done in within the function
        let lookupOperation = LookupOperationWrapper(minInFxP, outBits, numSwapBits, _, _);


        return FunctionWithLookupTable(
            pIn, qIn,
            pOut, qOut,
            lookupOperation
        );
    }

    /// # Summary
    /// Creates a lookup table operation. This operation will require the minimum input value as a FixedPoint register, 
    /// the list of output values in bits,the FixedPoint register with the input value and the FixedPoint register that 
    /// will store the output value. Note that this imples that the bit size requirement of these registers are pre-computed 
    /// beforehand
    ///
    /// # Input
    /// ## minInFxp
    /// The minimum possible value of the input to the lookup table
    /// ## outBits
    /// The list of output values in bits in order, where the first value is the output for the smallest input value and
    /// the last value is the output for the largest input value
    /// ## numSwapBits
    /// The number of bits of the input register that will be used in the SWAP section of the circuits. Another way of looking
    /// at this is that in step in the SELECT section of the circuit in Fig 1c of arXiv:1812.00954, we will encode 2^numSawpBits
    /// encoded
    /// ## input
    /// Qubit FixedPoint register containing input values
    /// ## output
    /// Qubit FixedPoint register containing where output values will be stored
    internal operation LookupOperationWrapper(minInFxP: Double, outBits: Bool[][], numSwapBits : Int, input: FixedPoint, output: FixedPoint) : Unit is Adj {

        let integerBitsIn = input::IntegerBits;
        let registerIn = input::Register;
        let fractionalBitsIn = Length(registerIn) - integerBitsIn;

        // We are now creating the lookup table. If the smallest value (i.e. minInFxP) happens to be 0, then we can just use
        // the Select operation which implements the lookup table in ##. However, if the minimum value is not 0, then we want to first subtract
        // it, because the lookup table always assumes that the miminum value is 00...0 and the maximum value is 11...1 in incrementing order,
        // so we are re-defining the minimum number as represented by 00...0 and hence subracting the minimum from our value.
        // (We add the minimum back after making the lookup table)
        within { // Currently we always uncompute the lookup table garbage qubits, but we can think of making an option to remove the uncomputation (and keep the garbage qubits)
            if minInFxP != 0.0 {
                // Make a new fixed point register to store the minimum vlaue
                use minRegister = Qubit[Length(registerIn)];
                let minInReg = FixedPoint(integerBitsIn, minRegister); //
                within {
                    PrepareFxP(minInFxP, minInReg); // Store minimum value in prepared register (automatically creates closest FxP approximation)
                } apply {
                    SubtractFxP(input, minInReg); // SubtractFxP(a, b) : a <- a - b
                }
            }
        } apply {
            SelectSwap(numSwapBits, outBits, input::Register, output::Register);
        }
    }

    /// # Summary
    /// Basically applies the ApplyPauliFromBitString function but with the extra check that the length of bits is the same
    /// as the number of qubits (so nothing can be implicitly ignored)
    internal operation WriteBits(bits: Bool[], qubitArray: Qubit[]): Unit is Adj + Ctl {
        EqualityFactI(Length(bits), Length(qubitArray), "Dimensions of bits and qubitArray should be the same");
        ApplyPauliFromBitString(PauliX, true, bits, qubitArray);
    }

    /// # Summary
    /// Helper function that creates an operator that takes in just 1 binary value input (i.e. a list 
    /// of booleans) and makes the circuit to apply paulis to create that binary value. We do this 
    /// so that we can use it as part of the Mapped function to be able to make a list of unitaries 
    /// given a list of binary numbers
    internal function MakeWriteBitsUnitary(bits : Bool[]) : Qubit[] => Unit is Adj + Ctl {
        return WriteBits(bits, _);
        
    }

    /// # Summary
    /// This opration makes the lookup table by using the multiplex unitary operator - the operator that implements
    /// different unitaries based on the value of the controlled bits. We just define each unitary as the set of
    /// PauliX gates that will make the output qubit correspond to the data bits.
    internal operation Select(data : Bool[][], addressRegister: Qubit[], outputRegister: Qubit[]) : Unit is Adj {

        let unitaries = Mapped(MakeWriteBitsUnitary, data);
        MultiplexOperations(unitaries, LittleEndian(addressRegister), outputRegister);
    }

    /// # Summary
    /// This operation makes the swap circuit. The outputRegisters are 2^l qubit registers, each of size m
    internal operation SwapDataOutputs(addressRegister: Qubit[], outputRegisters: Qubit[][]) : Unit is Adj {
        let l = Length(addressRegister);
        // For each input qubit we are using for swap qubits, we want to implement all the swaps
        for i in 0.. Length(addressRegister)-1{
            // for the ith input qubit, we need to have to swap qubit registers that are 2^i places apart, and we need to do this for every qubit that is 2^(i+1) qubits apart,
            // and we need to swap atotal of 2^l/2^(i+1) registers, where l is the number of output registers. E.g.
            // i=0 => swaps between registers (0,1), (2,3), (4,5),..., (2^l - 2, 2^l - 1)
            // i=1 => swaps between registers (0,2), (4,6), (8,10),..., (2^l - 4, 2^l - 2)
            // i=2 => swaps between registers (0,2^i), (2^(i+1), 2^(i+1) + 2^i),..., (2^l - 2^(i+1), 2^l - 2^(i+1) + 2^i)
            let innerStepSize = 2^i;  
            let outerStepSize = 2^(i+1);
            let numSwaps = 2^l/2^(i+1);
            use extraControls = Qubit[numSwaps-1];
            let fannedControls = [addressRegister[i]] + extraControls;
            within {
                ApplyToEachA(CNOT(addressRegister[i],_), extraControls); // Fanning out qubits to be able to do a fanned control
            } apply {
                for j in 0.. numSwaps-1{
                    ApplyMultiTargetSwap(fannedControls[j], outputRegisters[j*outerStepSize], outputRegisters[j*outerStepSize+innerStepSize]);
                }
            }
        }
    }
    
    /// # Summary   
    /// Implements an efficient control swap (with 2 CNOTs and one CCNOT)
    internal operation ApplyLowDepthCSWAP(control : Qubit, target1 : Qubit, target2 : Qubit) : Unit is Adj {
        use helper = Qubit();

        within {
            CNOT(target2, target1);
            ApplyLowDepthAnd(control, target1, helper); // this has T-depth 1, AND = CCNOT where target is in |0>
        } apply {
            CNOT(helper, target2);
        }
    }

    /// # Summary
    /// Implements a control swap two registers controlled on a single qubits. To be able to parallelize it, it will 
    /// fan out the control register and perform all the swaps in parallel
    internal operation ApplyMultiTargetSwap(control : Qubit, target1 : Qubit[], target2 : Qubit[]) : Unit is Adj {
        EqualityFactI(Length(target1), Length(target2), "The two qubit registers are of different sizes");

        use extraControls = Qubit[Length(target1) - 1];
        let fannedControls = [control] + extraControls;

        within {
            ApplyToEachA(CNOT(control, _), extraControls);
        } apply {
            for i in 0..Length(target1)-1{
                ApplyLowDepthCSWAP(fannedControls[i], target1[i], target2[i]);
            }
        }
    }

    /// # Summary
    /// Creates the select-swap circuit. Uses the most significant bits of `addressRegister` to use for SWAP network.
    /// Let n be `Length(addressRegister)`, the number of address bits, and let l be `numSwapBits`. Then the T-depth is 
    /// approximately O(2^{(n-l}) + l. If we want m output qubits, the number of additional ancilla qubits is m * 2^l:
    /// These due additional ancilla qubits come from the swap 
    ///
    /// # Input
    /// ## numSwapBits
    /// The number of bits of the input register that will be used in the SWAP section of the circuits. Another way of looking
    /// at this is that in step in the SELECT section of the circuit in Fig 1c of arXiv:1812.00954, we will encode 2^numSawpBits
    /// encoded.
    /// ## data
    /// The list of output values in bits in order, where the first value is the output for the value 00...0 and
    /// the last value is the output for the value the largest integer value
    /// ## addressRegister
    /// Input register of qubits
    /// ## outputRegister
    /// Output register where the output values will be stored    
    internal operation SelectSwap(numSwapBits : Int, data : Bool[][], addressRegister: Qubit[], outputRegister: Qubit[]) : Unit is Adj{
        
        let n = Length(addressRegister);

        // how many address bits do we need for `data`? We do this so that we optimize if Length(data) <= 2^(n-1) and don't have to use all the input qubits in the swap register
        let nRequired = Ceiling(Lg(IntAsDouble(Length(data))));
        Fact(nRequired <= n, "Too few address bits");
        let addressRegisterFitted = addressRegister[...nRequired - 1];

        // Probably numSwapBits == nRequired works, but has to be checked
        Fact(numSwapBits < nRequired, "Too many bits for SWAP network");

        if numSwapBits == 0 { // to not uncompute select if l=0
            Select(data, addressRegisterFitted, outputRegister);
        } else {
            let m = Length(outputRegister);
            // number of SWAP bits
            let l = numSwapBits;
            // number of SELECT bits
            let k = nRequired - numSwapBits;

            let addressRegisterParts = Partitioned([k, l], addressRegisterFitted); //creates two parts of the array - 1 for the select and 1 for the swaps

            use dataRegister = Qubit[m * 2^l]; //create one register with m*2^l qubits to have all the outputs stored

            // Now we want to create the data array for the select bit. This means that for the first array, we want all the outputs corresponding to the most sig bits = 00..00, then for the next array 00..01, then 00..10 etc. 
            let dataArray = CreatePaddedData(data, nRequired, m, k);

            // Now we split up the register in to chunks of m
            let chunkedDataRegister = Chunks(m, dataRegister); //Chunks(nElements, array): splits array into chunks, where each chunk has nElements, e.g. Chunks(2, [1, 2, 3, 4, 5, 6]) -> [[1, 2], [3, 4], [5, 6]]

            // Perform select swap with the (1) dataArray for the select outputs and (2) chunkedDataRegister for the swap targets. We can think about improving the efficiency of these two steps by thinking of the fanout inner and outer controls more carefully.
            within {
                Select(dataArray, addressRegisterParts[0], dataRegister); // apply select using first part of address registers
                SwapDataOutputs(addressRegisterParts[1], chunkedDataRegister); // apply swap using second part of address registers
            } apply {
                ApplyToEachA(CNOT, Zipped(chunkedDataRegister[0], outputRegister)); // apply CNOT from chunkedDataRegister[0] to outputRegister
            } 
        }
    }

    /// # Summary
    /// Helper function that creates a list of binary numbers were each number is a concatination of all the possible outputs of each
    /// select function (i.e. combining all the possible outputs into a single bitstring which will then be used in the swap bit of the
    /// select-swap)
    ///
    /// # Input
    /// ## data
    /// The list of output values in bits in order, where the first value is the output for the value 00...0 and
    /// the last value is the output for the value the largest integer value
    /// ## nRequired
    /// Minimum number of address bits required to ensure the largest data point can be stored
    /// ## m
    /// Size of final single output of lookup table
    /// ## k
    /// Number of bits from input that will be used in the select part
    internal function CreatePaddedData(data : Bool[][], nRequired : Int, m : Int, k : Int) : Bool[][] {
        let dataPadded = Padded(-2^nRequired, [false, size = m], data); // Padded so that we don't have problems for powers of 2
        mutable dataArray = [[], size = 2^k];
        for i in 0..2^k-1 {
            let range = RangeAsIntArray(i..2^k..2^nRequired-1); // RangeAsIntArray(range) takes a Range and returns Int[], e.g., 1..3 -> [1, 2, 3], 1..3..9 -> [1, 4, 7]. In our case, for i=0 -> [0, 2^k, 2*2^k...]
            set dataArray w/= i <- Flattened(Subarray(range, dataPadded)); // Subarray(indices, array) takes indices Int[] and returns subarray of array subject to indices
        } 
        return dataArray;
    }
}
