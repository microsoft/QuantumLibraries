// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Chemistry.JordanWigner {
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Extensions.Math;
    open Microsoft.Quantum.Chemistry;
    
    
    /// # Summary
    /// Format of data passed from C# to Q# to represent terms of the Hamiltonian.
    /// The meaning of the data represented is determined by the algorithm that receives it.
    newtype JWOptimizedHTerms = (HTerm[], HTerm[], HTerm[], HTerm[]);
    
    /// # Summary
    /// Format of data passed from C# to Q# to represent preparation of the initial state
    /// The meaning of the data represented is determined by the algorithm that receives it.
    newtype JordanWignerInputState = ((Double, Double), Int[]);
    
    /// # Summary
    /// Format of data passed from C# to Q# to represent all information for Hamiltonian simulation.
    /// The meaning of the data represented is determined by the algorithm that receives it.
    newtype JordanWignerEncodingData = (Int, JWOptimizedHTerms, JordanWignerInputState[], Double);
    
}


