// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {
    open Microsoft.Quantum.Extensions.Math;


    /// # Summary
    /// Represents a single primitive term in the set of all dynamical generators, e.g.
    /// Hermitian operators, for which there exists a map from that generator
    /// to time-evolution by that that generator, through "EvolutionSet". The first element
    /// (Int[], Double[]) is indexes that single term -- For instance, the Pauli string
    /// XXY with coefficient 0.5 would be indexed by ([1,1,2], [0.5]). Alternatively,
    /// Hamiltonians parameterized by a continuous variable, such as X cos φ + Y sin φ,
    /// might for instance be represented by ([], [φ]). The second
    /// element indexes the subsystem on which the generator acts on.
    ///
    /// # Remarks
    /// > [!WARNING]
    /// > The interpretation of an `GeneratorIndex` is not defined except
    /// > with reference to a particular set of generators.
    ///
    /// # Example
    /// Using  @"microsoft.quantum.canon.paulievolutionset", the operator
    /// $\pi X_2 X_5 Y_9$ is represented as:
    /// ```qsharp
    /// let index = GeneratorIndex(([1; 1; 2], [PI()]), [2; 5; 9]);
    /// ```
    ///
    /// # See Also
    /// - @"microsoft.quantum.canon.paulievolutionset"
    /// - @"microsoft.quantum.canon.fermionicevolutionset"
    newtype GeneratorIndex = ((Int[], Double[]), Int[]);

    // FIXME: unify this and the GateSet representation above with ContinousOracle.
    // FIXME: add an example using lookupfunction.
    /// # Summary
    /// From the view of a GeneratorSystem, a description of a Hamiltonian
    //  is a collection of GeneratorTerms. We iterate over this
    /// collection using a single-index integer, and the size of the
    /// collection is assumed to be known.
    ///
    /// # Remarks
    /// Instances of `GeneratorSystem` can be defined easily using the
    /// @"microsoft.quantum.canon.lookupfunction" function.
    ///
    /// # See Also
    /// - @"microsoft.quantum.canon.lookupfunction"
    newtype GeneratorSystem = (Int, (Int -> GeneratorIndex));

    /// # Summary
    /// Represents a time-dependent dynamical generator as a function
    /// from time to the value of the dynamical generator at that time.
    newtype GeneratorSystemTimeDependent = (Double -> GeneratorSystem);

    /// # Summary
    /// Returns a generator index consistent with the
    /// Hamiltonian $H = 0$ corresponding to the evolution $U(t) = \boldone$.
    ///
    /// # Input
    /// ## idxTerm
    /// This input is ignored, and is defined for consistency with the
    /// @"microsoft.quantum.canon.generatorsystem" user-defined type.
    ///
    /// # Output
    /// A generator index representing evolution under the Hamiltonian
    /// $H = 0$.
    function IdentityGeneratorIndex(idxTerm : Int) : GeneratorIndex {
        return GeneratorIndex(([0], [Float(0)]),[0]);
    }
    function IdentityGeneratorSystem() : GeneratorSystem {
        return  GeneratorSystem(0, IdentityGeneratorIndex);
    }
    function _IdentityGeneratorSystemTimeDependent(schedule: Double) : GeneratorSystem {
        return IdentityGeneratorSystem();
    }
    function IdentityGeneratorSystemTimeDependent() : GeneratorSystemTimeDependent {
        return GeneratorSystemTimeDependent(_IdentityGeneratorSystemTimeDependent);
    }

    /// Retrieves the number of terms in a GeneratorSystem
    function GetGeneratorSystemNTerms(generatorSystem : GeneratorSystem) : Int {
        let (nTerms, generatorIndexFunction) = generatorSystem;
        return nTerms;
    }
    /// Retrieves the number of terms in a GeneratorSystem
    function GetGeneratorSystemFunction(generatorSystem : GeneratorSystem) : (Int -> GeneratorIndex) {
        let (nTerms, generatorIndexFunction) = generatorSystem;
        return generatorIndexFunction;
    }
    
    // We should be able to do some algebra on representations of the system  
    /// Multiplies the coefficient index in a GeneratorIndex = ((Int[],Double[]), Int[])
    function MultiplyGeneratorIndex(multipler: Double, generatorIndex : GeneratorIndex) : GeneratorIndex{
        let ((idxTerms, idxDoubles), idxSystems) = generatorIndex;
        mutable idxDoublesOut = idxDoubles;
        set idxDoublesOut[0] = multipler * idxDoublesOut[0];
        return GeneratorIndex((idxTerms, idxDoublesOut), idxSystems);
    }
    /// Multiples the coefficient index of all terms in a GeneratorSystem
    function _MultiplyGeneratorSystem(multipler: Double, idxTerm: Int, generatorSystem : GeneratorSystem) : GeneratorIndex {
        let (nTerms, generatorIndexFunction) = generatorSystem;
        return MultiplyGeneratorIndex(multipler , generatorIndexFunction(idxTerm));
    }
    function MultiplyGeneratorSystem(multipler: Double, generatorSystem : GeneratorSystem) : GeneratorSystem {
        let nTerms =  GetGeneratorSystemNTerms(generatorSystem);
        return GeneratorSystem(nTerms, _MultiplyGeneratorSystem(multipler, _, generatorSystem));
    }

    /// Adds GeneratorSystems to create a new GeneratorSystem.
    function _AddGeneratorSystems(  idxTerm : Int, 
                                    nTermsA : Int, 
                                    nTermsB : Int, 
                                    generatorIndexFunctionA : (Int -> GeneratorIndex), 
                                    generatorIndexFunctionB : (Int -> GeneratorIndex)) : GeneratorIndex {
        if(idxTerm < nTermsA){
            return generatorIndexFunctionA(idxTerm);
        }
        else{
            return generatorIndexFunctionB(idxTerm - nTermsA);
        }
    }
    function AddGeneratorSystems(generatorSystemA: GeneratorSystem, generatorSystemB: GeneratorSystem) : GeneratorSystem {
        let nTermsA = GetGeneratorSystemNTerms(generatorSystemA);
        let nTermsB = GetGeneratorSystemNTerms(generatorSystemB);
        let generatorIndexFunctionA = GetGeneratorSystemFunction(generatorSystemA);
        let generatorIndexFunctionB = GetGeneratorSystemFunction(generatorSystemB);
        let generatorIndexFunction = _AddGeneratorSystems(_, nTermsA, nTermsB, generatorIndexFunctionA, generatorIndexFunctionB);
        return GeneratorSystem(nTermsA + nTermsB, generatorIndexFunction);
    }

    /// Create description of system that interpolates between two GeneratorSystems based on a schedule parameter in [0,1]
    function InterpolateGeneratorSystemsImpl(schedule: Double, generatorSystemStart: GeneratorSystem, generatorSystemEnd: GeneratorSystem) : GeneratorSystem
    {
        let sysStart = MultiplyGeneratorSystem((1.0 - schedule), generatorSystemStart);
        let sysEnd = MultiplyGeneratorSystem(schedule, generatorSystemEnd);
        return AddGeneratorSystems(sysStart, sysEnd);
    } 
    function InterpolateGeneratorSystems(generatorSystemStart: GeneratorSystem, generatorSystemEnd: GeneratorSystem) : GeneratorSystemTimeDependent
    {
        return GeneratorSystemTimeDependent(InterpolateGeneratorSystemsImpl(_, generatorSystemStart, generatorSystemEnd));
    } 

}
