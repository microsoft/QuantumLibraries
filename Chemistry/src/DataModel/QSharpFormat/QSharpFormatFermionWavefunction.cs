// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using Microsoft.Quantum.Simulation.Core;

using System;
using System.Linq;
using System.Collections.Generic;
using Microsoft.Quantum.Chemistry;
using Microsoft.Quantum.Chemistry.LadderOperators;
using Microsoft.Quantum.Chemistry.Broombridge;
using Microsoft.Quantum.Chemistry.Fermion;

namespace Microsoft.Quantum.Chemistry.QSharpFormat
{
    using Microsoft.Quantum.Chemistry.JordanWigner;
   
    /// <summary>
    /// Methods for converting electronic structure problem to data for consumption by Q#.
    /// </summary>
    public static partial class Convert
    {        
        /// <summary>
        /// Translate initial state to a format consumable by Q#.
        /// </summary>
        /// <param name="inputState">Initial state</param>
        /// <returns>Initial state in Q# format.</returns>
        public static (Int64, QArray<JordanWignerInputState>) ToQSharpFormat(this FermionWavefunction<int> inputState)
        {
            //return (0, new QArray<JordanWignerInputState>());
            
            if(inputState.Method == StateType.SparseMultiConfigurational)
            {
                return ((int)inputState.Method, InitialStateSparseMultiConfigural(inputState.MCFData));
            }
            else if(inputState.Method == StateType.UnitaryCoupledCluster)
            {
                return ((int)inputState.Method, InitialStateUnitaryCoupledCluster(inputState.UCCData));
            }
            else
            {
                throw new ArgumentException($"Selected quantum state {inputState.Method.ToString("G")} not recognized.");
                }           
        }

        /// <summary>
        /// Translate initial state specified by a superposition of <see cref="FermionTerm"/> acting on 
        /// the vacuum state to a format consumable by Q#
        /// </summary>
        /// /// <param name="term">List of sequences of creation operations acting on vacuum state and their coefficients.</param>
        /// <returns>Q# description of initial state</returns>
        public static QArray<JordanWignerInputState> InitialStateSparseMultiConfigural(
            SparseMultiCFWavefunction<int> wavefunction)
        {
            return new QArray<JordanWignerInputState>(
                wavefunction.Excitations.Select(o =>  InitialStatePrep(o.Value, o.Key, checkAnnihilation: true))
                // Todo add reference wavefunction when that functionality is implemented.
                );
        }

        /// <summary>
        /// Translate initial state specified by a <see cref="FermionTerm"/> acting on 
        /// the vacuum state to a format consumable by Q#
        /// </summary>
        /// <param name="term">Sequence of creation operations acting on vacuum state.</param>
        /// <param name="complexCoeff">coefficient of term.</param>
        /// <param name="checkAnnihilation">If true, throws an exception if annihilation operators act on the vacuum state.</param>
        /// <returns>Q# description of initial state</returns>
        public static JordanWignerInputState InitialStatePrep(System.Numerics.Complex complexCoeff, IndexOrderedSequence<int> term, bool checkAnnihilation)
        {
            // Place input in canonical order.
            var termCanonical = term;

            if (checkAnnihilation)
            {

                // There must be no annihilation operators the presence of
                // any of the in normal order destroys the vacuum state.
                if (termCanonical.Sequence.Where(o => o.Type == RaisingLowering.d).Count() > 0)
                {
                    throw new System.ArgumentException("Initial state cannot contain annihilation operators acting on the vacuum state.");
                }

            }

            // There must be no duplicate creation operators.
            if(termCanonical.Sequence.Distinct().Count() != termCanonical.Sequence.Count())
            {
                throw new System.ArgumentException("Initial state cannot contain duplicate creation operators with the same index.");
            }

            // Rescale input complex coefficient by fermion term coefficient,
            // including any sign changes.
            (double, double) complexOut = (complexCoeff.Real, complexCoeff.Imaginary);
            QArray<Int64> indices = new QArray<Int64>(term.Sequence.Select(o => (Int64)o.Index));
            var state = new JordanWignerInputState((complexOut, indices));
            return state;
        }


        /// <summary>
        /// Translate initial state specified by a unitary coupled-cluster operator acting on a single-reference state
        /// to a format consumable by Q#.
        /// </summary>
        /// <param name="terms">single-reference state and cluster operator terms.</param>
        /// <returns>Q# description of unitary coupled-cluster state.</returns>
        internal static QArray<JordanWignerInputState> InitialStateUnitaryCoupledCluster(
            UnitaryCCWavefunction<int> wavefunction)
        {
            // The last term is the reference state.
            var referenceState = wavefunction.Reference;
            var clusterOperator = wavefunction.Excitations;

            var stateElements = clusterOperator.Select(o => InitialStatePrep(o.Value, o.Key, checkAnnihilation: false)).ToList();
            stateElements.Add(InitialStatePrep(new System.Numerics.Complex(1.0, 0.0), referenceState, checkAnnihilation: true));

            var stateQSharp = new QArray<JordanWignerInputState>(stateElements);
            
            return stateQSharp;
        }
    }

}