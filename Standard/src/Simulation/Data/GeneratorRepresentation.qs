// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Simulation {
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Arrays;


    /// # Summary
    /// Represents a single primitive term in the set of all dynamical generators, e.g.
    /// Hermitian operators, for which there exists a map from that generator
    /// to time-evolution by that generator, through `EvolutionSet`.
    ///
    /// The first element
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
    /// Using <xref:Microsoft.Quantum.Simulation.PauliEvolutionSet>, the operator
    /// $\pi X_2 X_5 Y_9$ is represented as:
    /// ```qsharp
    /// let index = GeneratorIndex(([1, 1, 2], [PI()]), [2, 5, 9]);
    /// ```
    ///
    /// # See Also
    /// - Microsoft.Quantum.Simulation.EvolutionSet
    /// - Microsoft.Quantum.Simulation.PauliEvolutionSet
    newtype GeneratorIndex = ((Int[], Double[]), Int[]);

    /// # Summary
    /// Represents a collection of `GeneratorIndex`es.
    ///
    /// We iterate over this
    /// collection using a single-index integer, and the size of the
    /// collection is assumed to be known.
    ///
    /// # Remarks
    /// Instances of `GeneratorSystem` can be defined easily using the
    /// <xref:Microsoft.Quantum.Arrays.LookupFunction> function.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Arrays.LookupFunction
    newtype GeneratorSystem = (Int, (Int -> GeneratorIndex));

    /// # Summary
    /// Represents a time-dependent dynamical generator as a function
    /// from time to the value of the dynamical generator at that time.
    newtype TimeDependentGeneratorSystem = (Double -> GeneratorSystem);

    /// # Summary
    /// Represents evolution under a unitary operator.
    newtype Unitary = (Qubit[] => Unit is Adj + Ctl);

    /// # Summary
    /// Represents a unitary time-evolution operator.
    ///
    /// The first parameter is
    /// is duration of time-evolution, and the second parameter is the qubit
    /// register acted upon by the unitary.
    newtype EvolutionUnitary = ((Double, Qubit[]) => Unit is Adj + Ctl);

    /// # Summary
    /// Represents a set of gates that can be readily implemented and used
    /// to implement simulation algorithms.
    ///
    /// Elements in the set are indexed
    /// by a <xref:Microsoft.Quantum.Simulation.GeneratorIndex>,
    /// and each set is described by a function
    /// from `GeneratorIndex` to <xref:Microsoft.Quantum.Simulation.EvolutionUnitary>,
    /// which are operations
    /// parameterized by a real number representing time
    newtype EvolutionSet = (GeneratorIndex -> EvolutionUnitary);

    /// # Summary
    /// Represents a dynamical generator as a set of simulatable gates and
    /// an expansion in terms of that basis.
    ///
    /// Last parameter for number of terms.
    newtype EvolutionGenerator = (EvolutionSet, GeneratorSystem);

    /// # Summary
    /// Represents a time-dependent dynamical generator.
    ///
    /// The `Double`
    /// parameter is a schedule in $[0, 1]$.
    newtype EvolutionSchedule = (EvolutionSet, (Double -> GeneratorSystem));


    /// # Summary
    /// Returns a generator index consistent with the zero
    /// Hamiltonian, `H = 0`, which corresponds to the identity evolution operation.
    ///
    /// # Input
    /// ## idxTerm
    /// This input is ignored, and is defined for consistency with the
    /// <xref:Microsoft.Quantum.Simulation.GeneratorSystem> user-defined type.
    ///
    /// # Output
    /// A generator index representing evolution under the Hamiltonian
    /// $H = 0$.
    function IdentityGeneratorIndex (idxTerm : Int) : GeneratorIndex
    {
        return GeneratorIndex(([0], [0.0]), [0]);
    }


    /// # Summary
    /// Returns a generator system consistent with the zero
    /// Hamiltonian `H = 0`, which corresponds to the identity evolution operation.
    ///
    /// # Output
    /// A generator system representing evolution under the Hamiltonian
    /// $H = 0$.
    function IdentityGeneratorSystem () : GeneratorSystem {
        return GeneratorSystem(0, IdentityGeneratorIndex);
    }


    /// # Summary
    /// Returns a generator system consistent with the
    /// Hamiltonian `H(s) = 0`, where `s` is a schedule parameter.
    ///
    /// # Input
    /// ## schedule
    /// This input is ignored, and is defined for consistency with the
    /// <xref:Microsoft.Quantum.Simulation.TimeDependentGeneratorSystem> user-defined type.
    ///
    /// # Output
    /// A generator system representing evolution under the Hamiltonian
    /// $H(s) = 0$ for all $s$.
    internal function _IdentityTimeDependentGeneratorSystem (schedule : Double) : GeneratorSystem
    {
        return IdentityGeneratorSystem();
    }


    /// # Summary
    /// Returns a time-dependent generator system consistent with the
    /// Hamiltonian `H(s) = 0`.
    ///
    /// # Output
    /// A time dependent generator system representing evolution under the Hamiltonian
    /// $H(s) = 0$ for all $s$.
    function IdentityTimeDependentGeneratorSystem () : TimeDependentGeneratorSystem
    {
        return TimeDependentGeneratorSystem(_IdentityTimeDependentGeneratorSystem);
    }


    /// # Summary
    /// Retrieves the number of terms in a `GeneratorSystem`.
    ///
    /// # Input
    /// ## generatorSystem
    /// The `GeneratorSystem` of interest.
    ///
    /// # Output
    /// The number of terms in a `GeneratorSystem`.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Simulation.GeneratorSystem
    function GetGeneratorSystemNTerms (generatorSystem : GeneratorSystem) : Int {
        let (nTerms, generatorIndexFunction) = generatorSystem!;
        return nTerms;
    }


    /// # Summary
    /// Retrieves the `GeneratorIndex` function in a `GeneratorSystem`.
    ///
    /// # Input
    /// ## generatorSystem
    /// The `GeneratorSystem` of interest.
    ///
    /// # Output
    /// An function that indexes each `GeneratorIndex` term in a Hamiltonian.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Simulation.GeneratorIndex
    /// - Microsoft.Quantum.Simulation.GeneratorSystem
    function GetGeneratorSystemFunction (generatorSystem : GeneratorSystem) : (Int -> GeneratorIndex)
    {
        let (nTerms, generatorIndexFunction) = generatorSystem!;
        return generatorIndexFunction;
    }


    // We should be able to do some algebra on representations of the system

    /// # Summary
    /// Multiplies the coefficient in a `GeneratorIndex`.
    ///
    /// # Input
    /// ## multiplier
    /// The multiplier on the coefficient.
    /// ## generatorIndex
    /// The `GeneratorIndex` to be multiplied.
    ///
    /// # Output
    /// A `GeneratorIndex` representing a term with coefficient a factor
    /// `multiplier` larger.
    ///
    /// # Example
    /// ```qsharp
    /// let gen = GeneratorIndex(([1,2,3],[coeff]),[1,2,3]);
    /// let ((idxPaulis, idxDoubles), idxQubits) = MultiplyGeneratorIndex(multiplier, gen);
    /// // idxDoubles[0] == multiplier * coeff;
    /// ```
    ///
    /// # See Also
    /// - Microsoft.Quantum.Simulation.GeneratorIndex
    function MultiplyGeneratorIndex (multiplier : Double, generatorIndex : GeneratorIndex)
    : GeneratorIndex {
        let ((idxTerms, idxDoubles), idxSystems) = generatorIndex!;
        mutable idxDoublesOut = idxDoubles;
        set idxDoublesOut w/= 0 <- multiplier * idxDoublesOut[0];
        return GeneratorIndex((idxTerms, idxDoublesOut), idxSystems);
    }


    /// # Summary
    /// Multiplies the coefficient of all terms in a `GeneratorSystem`.
    ///
    /// # Remarks
    /// This is an intermediate step and should not be called.
    internal function _MultiplyGeneratorSystem (multiplier : Double, idxTerm : Int, generatorSystem : GeneratorSystem)
    : GeneratorIndex {
        let (nTerms, generatorIndexFunction) = generatorSystem!;
        return MultiplyGeneratorIndex(multiplier, generatorIndexFunction(idxTerm));
    }


    /// # Summary
    /// Multiplies the coefficient of all terms in a `GeneratorSystem`.
    ///
    /// # Input
    /// ## multiplier
    /// The multiplier on the coefficient.
    /// ## generatorSystem
    /// The `GeneratorSystem` to be multiplied.
    ///
    /// # Output
    /// A `GeneratorSystem` representing a system with coefficients a factor
    /// `multiplier` larger.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Simulation.GeneratorSystem
    function MultiplyGeneratorSystem (multiplier : Double, generatorSystem : GeneratorSystem) : GeneratorSystem
    {
        let nTerms = GetGeneratorSystemNTerms(generatorSystem);
        return GeneratorSystem(nTerms, _MultiplyGeneratorSystem(multiplier, _, generatorSystem));
    }


    /// # Summary
    /// Adds two `GeneratorSystem`s to create a new `GeneratorSystem`.
    ///
    /// # Remarks
    /// This is an intermediate step and should not be called.
    internal function _AddGeneratorSystems (idxTerm : Int, nTermsA : Int, nTermsB : Int, generatorIndexFunctionA : (Int -> GeneratorIndex), generatorIndexFunctionB : (Int -> GeneratorIndex)) : GeneratorIndex
    {
        if (idxTerm < nTermsA)
        {
            return generatorIndexFunctionA(idxTerm);
        }
        else
        {
            return generatorIndexFunctionB(idxTerm - nTermsA);
        }
    }


    /// # Summary
    /// Adds two `GeneratorSystem`s to create a new `GeneratorSystem`.
    ///
    /// # Input
    /// ## generatorSystemA
    /// The first `GeneratorSystem`.
    /// ## generatorSystemB
    /// The second `GeneratorSystem`.
    ///
    /// # Output
    /// A `GeneratorSystem` representing a system that is the sum of the
    /// input generator systems.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Simulation.GeneratorSystem
    function AddGeneratorSystems (generatorSystemA : GeneratorSystem, generatorSystemB : GeneratorSystem) : GeneratorSystem
    {
        let nTermsA = GetGeneratorSystemNTerms(generatorSystemA);
        let nTermsB = GetGeneratorSystemNTerms(generatorSystemB);
        let generatorIndexFunctionA = GetGeneratorSystemFunction(generatorSystemA);
        let generatorIndexFunctionB = GetGeneratorSystemFunction(generatorSystemB);
        let generatorIndexFunction = _AddGeneratorSystems(_, nTermsA, nTermsB, generatorIndexFunctionA, generatorIndexFunctionB);
        return GeneratorSystem(nTermsA + nTermsB, generatorIndexFunction);
    }


    // Create description of system that interpolates between two GeneratorSystems based on a schedule parameter in [0,1]
    //
    /// # Summary
    /// Linearly interpolates between two `GeneratorSystems` according to a
    /// schedule parameter `s` between 0 and 1 (inclusive).
    ///
    /// # Input
    /// ## schedule
    /// A schedule parameter $s\in[0,1]$.
    /// ##  generatorSystemStart
    /// The start `GeneratorSystem`.
    /// ## generatorSystemEnd
    /// The end `GeneratorSystem`.
    ///
    /// # Output
    /// A `GeneratorSystem` representing a system that is the sum of the
    /// input generator systems, with weight $(1-s)$ on `generatorSystemStart`
    /// and weight $s$ on `generatorSystemEnd`.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Simulation.GeneratorSystem
    internal function InterpolateGeneratorSystemsImpl (schedule : Double, generatorSystemStart : GeneratorSystem, generatorSystemEnd : GeneratorSystem) : GeneratorSystem
    {
        let sysStart = MultiplyGeneratorSystem(1.0 - schedule, generatorSystemStart);
        let sysEnd = MultiplyGeneratorSystem(schedule, generatorSystemEnd);
        return AddGeneratorSystems(sysStart, sysEnd);
    }


    /// # Summary
    /// Returns a `TimeDependentGeneratorSystem` representing the linear
    /// interpolation between two `GeneratorSystem`s.
    ///
    /// # Input
    /// ##  generatorSystemStart
    /// The start `GeneratorSystem`.
    /// ## generatorSystemEnd
    /// The end `GeneratorSystem`.
    ///
    /// # Output
    /// A `TimeDependentGeneratorSystem` representing a system that is the
    /// sum of the input generator systems, with weight $(1-s)$ on
    /// `generatorSystemStart` and weight $s$ on `generatorSystemEnd`.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Simulation.GeneratorSystem
    /// - Microsoft.Quantum.Simulation.TimeDependentGeneratorSystem
    function InterpolateGeneratorSystems (generatorSystemStart : GeneratorSystem, generatorSystemEnd : GeneratorSystem) : TimeDependentGeneratorSystem
    {
        return TimeDependentGeneratorSystem(InterpolateGeneratorSystemsImpl(_, generatorSystemStart, generatorSystemEnd));
    }

    /// # Summary
    /// Adds multiple `GeneratorSystem`s to create a new GeneratorSystem.
    ///
    /// # Input
    /// ## generatorSystems
    /// An array of type `GeneratorSystem[]`.
    ///
    /// # Output
    /// A `GeneratorSystem` representing a system that is the sum of the
    /// input generator systems.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Simulation.GeneratorSystem
    function SumGeneratorSystems(generatorSystems: GeneratorSystem[]) : GeneratorSystem {
        return Fold(AddGeneratorSystems, IdentityGeneratorSystem(), generatorSystems);
    }

}
