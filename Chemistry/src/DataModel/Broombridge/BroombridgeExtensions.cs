// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using System;
using System.Linq;
using System.Collections.Generic;

using Microsoft.Quantum.Chemistry.OrbitalIntegrals;
using Microsoft.Quantum.Chemistry.Fermion;
using Microsoft.Quantum.Chemistry.LadderOperators;
using Microsoft.Quantum.Chemistry.Generic;

namespace Microsoft.Quantum.Chemistry.Broombridge
{
    /// <summary>
    /// Extensions for converting orbital integrals to fermion terms.
    /// </summary>
    public static partial class Extensions
    {
        /// <summary>
        /// Builds Hamiltonian from Broombridge if data is available.
        /// </summary>
        public static OrbitalIntegralHamiltonian ToOrbitalIntegralHamiltonian(
            this CurrentVersion.ProblemDescription broombridge)
        {

            // Add the identity terms
            var identityterm = broombridge.CoulombRepulsion.Value + broombridge.EnergyOffset.Value;
            var hamiltonian = ToOrbitalIntegralHamiltonian(broombridge.Hamiltonian);
            hamiltonian.Add(new OrbitalIntegral(), identityterm);
            return hamiltonian;
        }


        /// <summary>
        /// Builds Hamiltonian from Broombridge orbital integral data.
        /// </summary>
        internal static OrbitalIntegralHamiltonian ToOrbitalIntegralHamiltonian(
            this DataStructures.DataStructure.HamiltonianData hamiltonianData)
        {
            var hamiltonian = new OrbitalIntegralHamiltonian();

            // This will convert from Broombridge 1-indexing to 0-indexing.
            hamiltonian.Add
                (hamiltonianData.OneElectronIntegrals.Values
                .Select(o => new OrbitalIntegral(o.Item1
                .Select(k => (int)(k - 1)), o.Item2, OrbitalIntegral.Convention.Mulliken)
                .ToCanonicalForm())
                .Distinct());

            // This will convert from Broombridge 1-indexing to 0-indexing.
            // This will convert to Dirac-indexing.
            hamiltonian.Add
                (hamiltonianData.TwoElectronIntegrals.Values
                .Select(o => new OrbitalIntegral(o.Item1
                .Select(k => (int)(k - 1)), o.Item2, OrbitalIntegral.Convention.Mulliken)
                .ToCanonicalForm())
                .Distinct());
            
            return hamiltonian;
        }


        /*
        internal static WavefunctionFermionMCF ToWavefunctionFermionMCF(
            this CurrentVersion.State inputState,
            SpinOrbital.IndexConvention indexConvention,
            int nOrbitals
            )
        {

            return ToWavefunctionFermionMCF(Deserializers.ParseInputState(inputState.Superposition), indexConvention, nOrbitals);
            
                
        }*/
        /*
        internal static WavefunctionFermionMCF ToWavefunctionFermionMCF(
            ((double, double), (int, Spin, RaisingLowering)[])[] inputState,
            SpinOrbital.IndexConvention indexConvention,
            int nOrbitals
            )
        {
            Func<(int, Spin, RaisingLowering), LadderOperator> ToLadderOperator = x
                => new LadderOperator(x.Item3, new SpinOrbital(x.Item1, x.Item2).ToInt(indexConvention, nOrbitals));

            Func<(int, Spin, RaisingLowering)[], WavefunctionFermionSCF> ToWavefunctionFermionSCF = x
                => new WavefunctionFermionSCF(new LadderSequence(x.Select(o => ToLadderOperator(o))));

            var state = inputState
                .Select(o => (o.Item1, ToWavefunctionFermionSCF(o.Item2))).ToList();

            return new WavefunctionFermionMCF()
            {
                Superposition = state,
                reference = new WavefunctionFermionSCF()
            };

        }
        */

        // This reads the number of orbitals from Broombridge
        public static Dictionary<string, InputState> ToWavefunctions(
            this CurrentVersion.ProblemDescription problemDescription,
            SpinOrbital.IndexConvention indexConvention)
        {
            var nOrbitals = problemDescription.NOrbitals;
            var states = problemDescription.InitialStates.Select(o => o.ToWavefunction(indexConvention, nOrbitals)).ToDictionary(o => o.Label, o => o);
            return states;
        }

        internal static ((double,double),IndexOrderedLadderSequence) ToSingleState(
            ((double, double), LadderSequence) superposition
            )
        {
            var ((amplitude, im), sequence) = superposition;
            var output = new IndexOrderedLadderSequence();
            var states = sequence.ToIndexOrder();
            var noAnnihilationTerm = states.Where(o => !(o.Sequence.Select(x => x.Type)).Contains(RaisingLowering.d));

            if(noAnnihilationTerm.Count() == 1)
            {
                output = noAnnihilationTerm.Single();
                amplitude *= (double)output.Coefficient;
            }
            return ((amplitude, 0.0), output);
        }

        internal static List<((double, double), IndexOrderedLadderSequence)> ToSingleState(
            ((double, double), (int, Spin, RaisingLowering)[])[] superposition,
            SpinOrbital.IndexConvention indexConvention,
            int nOrbitals
            )
        {
            Func<(int, Spin, RaisingLowering), LadderOperator> ToLadderOperator = x
                 => new LadderOperator(x.Item3, new SpinOrbital(x.Item1, x.Item2).ToInt(indexConvention, nOrbitals));

            Func<(int, Spin, RaisingLowering)[], LadderSequence> MakeState = x
                => new LadderSequence(x.Select(o => ToLadderOperator(o)));

            var output = superposition
                .Select(o => ToSingleState((o.Item1, MakeState(o.Item2))))
                .ToList();

            return output;
        }


        // Todo: Deserializer should refer spin orbital indices.
        internal static InputState ToWavefunction(
            this CurrentVersion.State initialState, 
            SpinOrbital.IndexConvention indexConvention,
            int nOrbitals)
        {


            Func<(int, Spin, RaisingLowering), LadderOperator> ToLadderOperator = x
                => new LadderOperator(x.Item3, new SpinOrbital(x.Item1, x.Item2).ToInt(indexConvention, nOrbitals));

            Func<(int, Spin, RaisingLowering)[], IndexOrderedLadderSequence> ToIndexOrder = x
                => new LadderSequence(x.Select(o => ToLadderOperator(o))).ToIndexOrder().First();

            var state = new InputState();
            state.TypeOfState = Deserializers.ParseInitialStateMethod(initialState.Method);
            state.Label = initialState.Label;
            if (state.TypeOfState == StateType.SparseMultiConfigurational)
            {
                ((double, double), (int, Spin, RaisingLowering)[])[]  deserializedState = Deserializers.ParseInputState(initialState.Superposition);

                state.Superposition = ToSingleState(deserializedState, indexConvention, nOrbitals);
            }
            else if (state.TypeOfState == StateType.UnitaryCoupledCluster)
            {
                var referenceState = Deserializers.ParseInputState(initialState.ClusterOperator.Reference);
                var reference = (referenceState.Item1, ToIndexOrder(referenceState.Item2));

                var oneBodyTerms = initialState.ClusterOperator.OneBodyAmplitudes.Select(o => Deserializers.ParseUnitaryCoupledClusterInputState(o))
                    .Select(o => (o.Item1, ToIndexOrder(o.Item2))).ToList();
                var twoBodyTerms = initialState.ClusterOperator.TwoBodyAmplitudes.Select(o => Deserializers.ParseUnitaryCoupledClusterInputState(o))
                    .Select(o => (o.Item1, ToIndexOrder(o.Item2))).ToList();
                var clusterTerms = oneBodyTerms.Concat(twoBodyTerms).ToList();
                // The last term is the reference state.
                clusterTerms.Add(reference);
                state.Superposition = clusterTerms;
                
            }
            else
            {
                throw new System.ArgumentException($"initial state `{state.Label}` is not recognized or implemented.");
            }
            return state;
        }
    }
}



