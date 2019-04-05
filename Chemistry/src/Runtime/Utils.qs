// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Chemistry {
    open Microsoft.Quantum.Simulation;
    open Microsoft.Quantum.Primitive;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Extensions.Math;
    
    
    /// # Summary
    /// Format of data passed from C# to Q# to represent a term of the Hamiltonian.
    /// The meaning of the data represented is determined by the algorithm that receives it.
    newtype HTerm = (Int[], Double[]);
    
    
    /// # Summary
    /// Converts a Hamiltonian term in `HTerm` data format to a GeneratorIndex.
    ///
    /// # Input
    /// ## term
    /// Input data in `HTerm` format.
    /// ## termType
    /// Additional information added to GeneratorIndex.
    ///
    /// # Output
    /// A GeneratorIndex representing a Hamiltonian term represented by `term`,
    /// together with additional information added by `termType`.
    function HTermToGenIdx (term : HTerm, termType : Int[]) : GeneratorIndex {
        
        let (idxFermions, coeff) = term!;
        return GeneratorIndex((termType, coeff), idxFermions);
    }
    
    
    /// # Summary
    /// Converts an index to a Hamiltonian term in `HTerm[]` data format to a GeneratorIndex.
    ///
    /// # Input
    /// ## data
    /// Input data in `HTerm[]` format.
    /// ## termType
    /// Additional information added to GeneratorIndex.
    /// ## idx
    /// Index to a term of the Hamiltonian
    ///
    /// # Output
    /// A GeneratorIndex representing a Hamiltonian term represented by `data[idx]`,
    /// together with additional information added by `termType`.
    function HTermsToGenIdx (data : HTerm[], termType : Int[], idx : Int) : GeneratorIndex {
        
        return HTermToGenIdx(data[idx], termType);
    }
    
    
    /// # Summary
    /// Converts a Hamiltonian in `HTerm[]` data format to a GeneratorSystem.
    ///
    /// # Input
    /// ## data
    /// Input data in `HTerm[]` format.
    /// ## termType
    /// Additional information added to GeneratorIndex.
    ///
    /// # Output
    /// A GeneratorSystem representing a Hamiltonian represented by the input `data`.
    function HTermsToGenSys (data : HTerm[], termType : Int[]) : GeneratorSystem {
        
        //return GeneratorSystem(Length(data), Compose(HTermToGenIdx(_, termType),LookupFunction(data)));
        return GeneratorSystem(Length(data), HTermsToGenIdx(data, termType, _));
    }
    
    
    /// # Summary
    /// Checks whether a `Double` number is not approximately zero.
    ///
    /// # Input
    /// ## number
    /// Number to be checked
    ///
    /// # Output
    /// Returns true if `number` has an absolute value greater than `1e-15`.
    function IsNotZero (number : Double) : Bool {
        
        if (AbsD(number) > PowD(10.0, -15.0)) {
            return true;
        }
        else {
            return false;
        }
    }
    
}


