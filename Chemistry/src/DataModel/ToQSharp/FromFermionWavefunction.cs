// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using Microsoft.Quantum.Simulation.Core;

using System;
using System.Linq;
using System.Collections.Generic;
using Microsoft.Quantum.Chemistry;

namespace Microsoft.Quantum.Chemistry
{
    using Microsoft.Quantum.Chemistry.JordanWigner;   

    public partial class JordanWignerEncoding
    {
        /// <summary>
        /// List of Hamiltonian input states.
        /// </summary>
        public Dictionary<string, FermionHamiltonian.InputState> InputStates;
        
        /// <summary>
        /// Translate initial state to a format consumable by Q#.
        /// </summary>
        /// <param name="inputState">Initial state</param>
        /// <returns>Initial state in Q# format.</returns>
        internal (int, QArray<JordanWignerInputState>) InputStateToQSharp(FermionHamiltonian.InputState inputState)
        {
            if(inputState.type == FermionHamiltonian.StateType.SparseMultiConfigurational)
            {
                return ((int)inputState.type, InitialStateSparseMultiConfigural(inputState.Superposition));
            }
            else if(inputState.type == FermionHamiltonian.StateType.UnitaryCoupledCluster)
            {
                return ((int)inputState.type, InitialStateUnitaryCoupledCluster(inputState.Superposition));
            }
            else
            {
                return InputStateToQSharp(InputStates["Greedy"]);
            }            
        }

        /// <summary>
        /// Translate initial state specified by a <see cref="FermionTerm"/> acting on 
        /// the vacuum state to a format consumable by Q#
        /// </summary>
        /// <param name="term">Sequence of creation operations acting on vacuum state.</param>
        /// <param name="complexCoeff"></param>
        /// <returns>Q# description of initial state</returns>
        public JordanWignerInputState InitialStatePrep((Double, Double) complexCoeff, FermionTerm term)
        {
            // Place input in canonical order.
            var termCanonical = term;
            termCanonical.ToSpinOrbitalCanonicalOrder();

            // There must be no annihilation operators the presence of
            // any of the in normal order destroys the vacuum state.
            if (termCanonical.CreationAnnihilationIndices.Contains(0))
            {
                throw new System.ArgumentException("Initial state cannot contain annihilation operators acting on the vacuum state.");
            }
            
            // There must be no duplicate creation operators.
            if(termCanonical.SpinOrbitalIndices.Distinct().Count() != termCanonical.SpinOrbitalIndices.Count())
            {
                throw new System.ArgumentException("Initial state cannot contain duplicate creation operators with the same index.");
            }

            // Rescale input complex coefficient by fermion term coefficient,
            // including any sign changes.
            (Double, Double) complexOut = (term.coeff * complexCoeff.Item1, complexCoeff.Item2 );
            var state = new JordanWignerInputState((complexOut, new QArray<Int64>(term.SpinOrbitalIndices.ToInts(NOrbitals))));
            return state;
        }

        /// <summary>
        /// Translate initial state specified by a superposition of <see cref="FermionTerm"/> acting on 
        /// the vacuum state to a format consumable by Q#
        /// </summary>
        /// /// <param name="term">List of sequences of creation operations acting on vacuum state and their coefficients.</param>
        /// <returns>Q# description of initial state</returns>
        public QArray<JordanWignerInputState> InitialStateSparseMultiConfigural(IEnumerable<((Double, Double) complexCoeff, FermionTerm term)> terms)
        {
            return new QArray<JordanWignerInputState>(terms.Select(o => InitialStatePrep(o.complexCoeff, o.term)));
        }

        /// <summary>
        /// Translate initial state specified by a unitary coupled-cluster operator acting on a single-reference state
        /// to a format consumable by Q#.
        /// </summary>
        /// <param name="terms">single-reference state and cluster operator terms.</param>
        /// <returns>Q# description of unitary coupled-cluster state.</returns>
        internal QArray<JordanWignerInputState> InitialStateUnitaryCoupledCluster(IEnumerable<((Double, Double) complexCoeff, FermionTerm term)> terms)
        {
            // The last term is the reference state.
            var referenceState = (terms.Last());
            var clusterOperator = terms.Take(terms.Count() - 1);
            var stateQSharp = new QArray<JordanWignerInputState>(clusterOperator.Select(o =>
            new JordanWignerInputState(
                (((Double) o.complexCoeff.Item1, (Double) o.complexCoeff.Item2),
                 new QArray<Int64>(o.term.SpinOrbitalIndices.ToInts(NOrbitals))))));
            stateQSharp.Add(InitialStatePrep(referenceState.complexCoeff, referenceState.term));
            return stateQSharp;
        }
    }

}