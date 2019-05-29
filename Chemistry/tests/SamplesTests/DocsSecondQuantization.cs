

// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

// This test ensures that any chemistry library syntax changes
// that affect samples are detected.

#region Using Statements
// We will need several different libraries in this sample.
// Here, we expose these libraries to our program using the
// C# "using" statement, similar to the Q# "open" statement.

// We will use the data model implemented by the Quantum Development Kit Chemistry
// Libraries. This model defines what a fermionic Hamiltonian is, and how to
// represent Hamiltonians on disk.
using Microsoft.Quantum.Chemistry.OrbitalIntegrals;
using Microsoft.Quantum.Chemistry.Fermion;
using Microsoft.Quantum.Chemistry.LadderOperators;


// The System namespace provides a number of useful built-in
// types and methods that we'll use throughout this sample.
using System;

// We use this for convnience functions for manipulation arrays.
using System.Linq;

//
using Xunit;
#endregion

namespace Microsoft.Quantum.Chemistry.Tests.Docs
{
    public static class SecondQuantization
    {
        [Fact]
        static void MakeSpinOrbital()
        {
            // First, we assign an orbital index, say `5`. Note that we use 0-indexing,
            // so this is the 6th orbital.
            var orbitalIdx = 5;

            // Second, we assign a spin index, say `Spin.u` for spin up or `Spin.d` for spin down.
            var spin = Spin.d;

            // The spin-orbital (5, ↓) is then
            var spinOrbital0 = new SpinOrbital(orbitalIdx, spin);

            // A tuple `(int, Spin)` is also implicitly recognized as a spin-orbital.
            (int, Spin) tuple = (orbitalIdx, spin);

            // We explicitly specify the type of `spinOrbital1` to demonstrate
            // the implicit cast to `SpinOrbital`.
            SpinOrbital spinOrbital1 = tuple;
        }

        [Fact]
        static void SpinOrbitalToInt()
        {
            // Let us use the spin orbital created in the previous snippet.
            var spinOrbital = new SpinOrbital(5, Spin.d);

            // Let us set the total number of orbitals to be say, `7`.
            var nOrbitals = 7;

            // This converts a spin-orbital index to a unique integer, in this case `12`,
            // using the formula `g(j,σ)`.
            var integerIndexHalfUp = spinOrbital.ToInt(IndexConvention.HalfUp, nOrbitals);

            // This converts a spin-orbital index to a unique integer, in this case `11`,
            // using the formula `h(j,σ)`.
            var integerIndexUpDown = spinOrbital.ToInt(IndexConvention.UpDown);

            // The default conversion uses the formula `h(j,σ)`, in this case `11`.
            var integerIndexDefault = spinOrbital.ToInt();

        }

        [Fact]
        static void LadderOperator()
        {
            // Let us use the spin orbital created in the previous snippet.
            var spinOrbitalInteger = new SpinOrbital(5, Spin.d).ToInt();

            // We specify either a creation or annihilation operator using 
            // the enumerable type `RaisingLowering.u` or `RaisingLowering.d`
            // respectively;
            var creationEnum = RaisingLowering.u;

            // The type representing a creation operator is then initialized 
            // as follows. Here, we index these operators with integers.
            // Hence we initialize the generic ladder operator with an
            // integer index type.
            var ladderOperator0 = new LadderOperator<int>(creationEnum, spinOrbitalInteger);

            // An alternate constructor for a LadderOperator instead uses
            // a tuple.
            var ladderOperator1 = new LadderOperator<int>((creationEnum, spinOrbitalInteger));
        }

        [Fact]
        static void FermionTerms()
        {
            // Let us initialize an array of tuples representing the
            // desired sequence of creation operators.
            var indices = new[] { (RaisingLowering.u, 1), (RaisingLowering.u, 2) };

            // We can convert this array of tuples to a sequence of ladder
            // operators using the `ToLadderSequence()` methods.
            var ladderSequences = indices.ToLadderSequence();

            // Sequences of ladder operators are quite general. For instance,
            // they could be bosonic operators, intead of fermionic operators.
            // We specialize them by constructing a `FermionTerm` representing 
            // a fermion creation operator on the index `2` followed by `1`.
            var fermionTerm = new FermionTerm(ladderSequences);
        }

        [Fact]
        static void NumberOperator()
        {
            // Let us use a new method to compactly create a sequence of ladder
            // operators. Note that we have ommitted specifying whether the 
            // operators are raising or lowering. In this case, the first half
            // will be raising operators, and the second half will be lowering 
            // operators.
            var indices = new[] { 1, 1 }.ToLadderSequence();

            // We now construct a `FermionTerm` representing an annihilation operator
            // on the index 1 followed by the creation operator on the index 1.
            var fermionTerm0 = new FermionTerm(indices);
        }

        [Fact]
        static void InequivalentLadderSequence()
        {
            // Let us initialize an array of tuples representing the
            // desired sequence of creation operators.
            var indices = new[] { (RaisingLowering.u, 1), (RaisingLowering.u, 2) };

            // We now construct a `LadderSequence` representing a creation operator
            // on the index 1 followed by 2, then a term with the reverse ordering.
            var laddderSeqeunce = indices.ToLadderSequence();
            var laddderSeqeunceReversed = indices.Reverse().ToLadderSequence();

            // The following Boolean is `false`.
            var equal = laddderSeqeunce == laddderSeqeunceReversed;
        }

        [Fact]
        static void EquivalentFermionTerm()
        {
            // We now construct two `FermionTerms` that are equivalent with respect to
            // anti-commutation up to a sign change.
            var fermionTerm0 = new FermionTerm(new[] { 0, 1, 1, 0 }.ToLadderSequence());
            var fermionTerm1 = new FermionTerm(new[] { 1, 0, 1, 0 }.ToLadderSequence());

            // The following Boolean is `true`.
            var sequenceEqual = fermionTerm0 == fermionTerm1;

            // The change in sign is not compared above, but is an internally tracked
            // property of `FermionTerm`.
            int sign0 = fermionTerm0.Coefficient;
            var sign1 = fermionTerm1.Coefficient;

            // The following Boolean is `false`.
            var signEqual = sign0 == sign1;
        }

        [Fact]
        static void EquivalentHermitianFermionTerm()
        {
            // We create a `FermionHamiltonian` object to store the fermion terms.
            var hamiltonian = new FermionHamiltonian();

            // We construct the terms to be added.
            var fermionTerm0 = new FermionTerm(new[] { 1, 0 }.ToLadderSequence());
            var fermionTerm1 = new FermionTerm(new[] { 0, 1 }.ToLadderSequence());

            // These fermion terms are not equal. The following Boolean is `false`.
            var sequenceEqual = fermionTerm0 == fermionTerm1;

            // However, these terms are equal under Hermitian symmetry.
            // We also take the opportunity to demonstrate equivalent constructors
            // for hermitian fermion terms
            var hermitianFermionTerm0 = new HermitianFermionTerm(fermionTerm0);
            var hermitianFermionTerm1 = new HermitianFermionTerm(new[] { 0, 1 });

            // These Hermitian fermion terms are equal. The following Boolean is `true`.
            var hermitianSequenceEqual = hermitianFermionTerm0 == hermitianFermionTerm1;

            // We add the terms to the Hamiltonian with the appropriate coefficient.
            // Note that these terms are identical.
            hamiltonian.Add(hermitianFermionTerm0, 1.0);
            hamiltonian.Add(hermitianFermionTerm1, 1.0);
        }

        [Fact]
        static void MakeHamiltonian()
        {
            // We create a `FermionHamiltonian` object to store the fermion terms.
            var hamiltonian = new FermionHamiltonian();

            // We construct the term to be added -- note the doubled coefficient.
            hamiltonian.Add(new HermitianFermionTerm(new[] { 1, 0 }), 2.0);
        }

        [Fact]
        static void MakeAnotherHamiltonian()
        {
            // We create a `FermionHamiltonian` object to store the fermion terms.
            var hamiltonian = new FermionHamiltonian();

            // We now create two Hermitian fermion terms that are equivalent with respect to
            // anti-commutation and Hermitian symmetry.
            var terms = new[] {
                (new[] { 0, 1, 1, 0 }, 1.0),
                (new[] { 1, 0, 1, 0 }, 1.0) }
            .Select(o => (new HermitianFermionTerm(o.Item1.ToLadderSequence()), o.Item2.ToDoubleCoeff()));

            // Now add `terms` to the Hamiltonian.
            hamiltonian.AddRange(terms);

            // There is only one unique term. `nTerms == 1` is `true`.
            var nTerms = hamiltonian.CountTerms();
        }
    }

    public static class Symmetries
    {
        [Fact]
        static void MakeOrbitalIntegral()
        {
            // We load the namespace containing orbital integral objects.
            // using Microsoft.Quantum.Chemistry.OrbitalIntegrals;

            // We create a `OrbitalIntegral` object to store a one-electron molecular 
            // orbital integral data.
            var oneElectronOrbitalIndices = new[] { 0, 1 };
            var oneElectronCoefficient = 1.0;
            var oneElectronIntegral = new OrbitalIntegral(oneElectronOrbitalIndices, oneElectronCoefficient);

            // This enumerates all one-electron integrals with the same coefficient --
            // an array of equivalent `OrbitalIntegral` instances is generated. In this
            //  case, there are two elements.
            var oneElectronIntegrals = oneElectronIntegral.EnumerateOrbitalSymmetries();

            // We create a `OrbitalIntegral` object to store a two-electron molecular orbital integral data.
            var twoElectronOrbitalIndices = new[] { 0, 1, 2, 3 };
            var twoElectronCoefficient = 0.123;
            var twoElectronIntegral = new OrbitalIntegral(twoElectronOrbitalIndices, twoElectronCoefficient);

            // This enumerates all two-electron integrals with the same coefficient -- 
            // an array of equivalent `OrbitalIntegral` instances is generated. In 
            // this case, there are 8 elements.
            var twoElectronIntegrals = twoElectronIntegral.EnumerateOrbitalSymmetries();

        }

        [Fact]
        static void MakeAnotherOrbitalIntegral()
        {
            // We create a `OrbitalIntegral` object to store a two-electron molecular
            //  orbital integral data.
            var twoElectronIntegral = new OrbitalIntegral(new[] { 0, 1, 2, 3 }, 0.123);

            // This enumerates all spin-orbital indices of the `FermionTerm`s in the 
            // Hamiltonian represented by this integral -- this is an array of array 
            // of `SpinOrbital` instances.
            var twoElectronSpinOrbitalIndices = twoElectronIntegral.EnumerateSpinOrbitals();

        }

        [Fact]
        static void MakeHamiltonian()
        {
            // We load the namespace containing fermion objects. This
            // example also uses LINQ queries.
            //using Microsoft.Quantum.Chemistry.Fermion;
            //using System.Linq;


            // We create a `OrbitalIntegral` object to store a two-electron molecular 
            // orbital integral data.
            var orbitalIntegral = new OrbitalIntegral(new[] { 0, 1, 2, 3 }, 0.123);

            // We create an `OrbitalIntegralHamiltonian` object to store the orbital integral
            // terms
            var orbitalIntegralHamiltonian = new OrbitalIntegralHamiltonian();
            orbitalIntegralHamiltonian.Add(orbitalIntegral);

            // We now convert the orbital integral representation to a fermion
            // representation. This also requires choosing a convention for 
            // mapping spin orbital indices to integer indices.
            var fermionHamiltonian = orbitalIntegralHamiltonian.ToFermionHamiltonian(IndexConvention.UpDown);

            // Alternatively, we we add orbital integrals directly to a fermion Hamiltonian
            // as follows. This automatically enumerates over all symmetries, and then
            // orders the `HermitianFermionTerm` instances in canonical order. We will need to
            // choose an indexing convention as well.
            fermionHamiltonian.AddRange(orbitalIntegral
                .ToHermitianFermionTerms(0, IndexConvention.UpDown)
                .Select(o => (o.Item1, o.Item2.ToDoubleCoeff())));
        }

    }
}