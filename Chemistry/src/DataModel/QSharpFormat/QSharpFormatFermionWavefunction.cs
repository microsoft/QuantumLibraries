// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using Microsoft.Quantum.Simulation.Core;

using System;
using System.Linq;
using System.Collections.Generic;
using Microsoft.Quantum.Chemistry;
using Microsoft.Quantum.Chemistry.Generic;
using Microsoft.Quantum.Chemistry.LadderOperators;

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
        public static (Int64, QArray<JordanWignerInputState>) ToQSharpFormat(this InputState inputState)
        {
            if(inputState.TypeOfState == StateType.SparseMultiConfigurational)
            {
                return ((int)inputState.TypeOfState, InitialStateSparseMultiConfigural(inputState.Superposition));
            }
            else if(inputState.TypeOfState == StateType.UnitaryCoupledCluster)
            {
                return ((int)inputState.TypeOfState, InitialStateUnitaryCoupledCluster(inputState.Superposition));
            }
            else
            {
                throw new ArgumentException($"Selected quantum state {inputState.TypeOfState.ToString("G")} not recognized.");
                return ((int)StateType.Default, null);
            }            
        }

        /// <summary>
        /// Translate initial state specified by a superposition of <see cref="FermionTerm"/> acting on 
        /// the vacuum state to a format consumable by Q#
        /// </summary>
        /// /// <param name="term">List of sequences of creation operations acting on vacuum state and their coefficients.</param>
        /// <returns>Q# description of initial state</returns>
        public static QArray<JordanWignerInputState> InitialStateSparseMultiConfigural(
            IEnumerable<((double, double) complexCoeff, IndexOrderedLadderSequence term)> terms)
        {
            return new QArray<JordanWignerInputState>(terms.Select(o => InitialStatePrep(o.complexCoeff, o.term)));
        }

        /// <summary>
        /// Translate initial state specified by a <see cref="FermionTerm"/> acting on 
        /// the vacuum state to a format consumable by Q#
        /// </summary>
        /// <param name="term">Sequence of creation operations acting on vacuum state.</param>
        /// <param name="complexCoeff"></param>
        /// <returns>Q# description of initial state</returns>
        public static JordanWignerInputState InitialStatePrep((double, double) complexCoeff, IndexOrderedLadderSequence term)
        {
            // Place input in canonical order.
            var termCanonical = term;

            // There must be no annihilation operators the presence of
            // any of the in normal order destroys the vacuum state.
            if (termCanonical.Sequence.Where(o => o.Type == RaisingLowering.d).Count() > 0)
            {
                throw new System.ArgumentException("Initial state cannot contain annihilation operators acting on the vacuum state.");
            }
            
            // There must be no duplicate creation operators.
            if(termCanonical.Sequence.Distinct().Count() != termCanonical.Sequence.Count())
            {
                throw new System.ArgumentException("Initial state cannot contain duplicate creation operators with the same index.");
            }

            // Rescale input complex coefficient by fermion term coefficient,
            // including any sign changes.
            (double, double) complexOut = ((double) term.Coefficient * complexCoeff.Item1, complexCoeff.Item2 );
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
            IEnumerable<((double, double) complexCoeff, IndexOrderedLadderSequence)> terms)
        {
            // The last term is the reference state.
            var referenceState = (terms.Last());
            var clusterOperator = terms.Take(terms.Count() - 1);

            var stateElements = clusterOperator.Select(o =>
            new JordanWignerInputState(
                (((Double)o.complexCoeff.Item1, (Double)o.complexCoeff.Item2),
                 new QArray<Int64>(o.Item2.Sequence.Select(y => (Int64)y.Index))))).ToList();
            stateElements.Add(InitialStatePrep(referenceState.complexCoeff, referenceState.Item2));

            var stateQSharp = new QArray<JordanWignerInputState>(stateElements);
            
            return stateQSharp;
        }
    }

}