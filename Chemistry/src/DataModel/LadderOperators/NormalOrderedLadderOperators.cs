// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using Microsoft.Quantum.Simulation.Core;

using System;
using System.Linq;
using System.Collections.Generic;
using Microsoft.Quantum.Chemistry;
using System.Numerics;

namespace Microsoft.Quantum.Chemistry
{
    
    public class NormalOrderedLadderOperators : LadderOperators
    {
        #region Constructors
        internal NormalOrderedLadderOperators() : base() { }

        /// <summary>
        /// Construct LadderOperators from another LadderOperators.
        /// </summary>
        /// <param name="setSequence">Sequence of ladder operators.</param>
        public NormalOrderedLadderOperators(NormalOrderedLadderOperators ladderOperators)
        {
            sequence = ladderOperators.sequence;
            coefficient = ladderOperators.coefficient;
        }

        /// <summary>
        /// Construct LadderOperators from sequence of ladder operators.
        /// </summary>
        /// <param name="setSequence">Sequence of ladder operators.</param>
        public NormalOrderedLadderOperators(LadderOperators ladderOperators) : base(ladderOperators)
        {
            ExceptionIfNotInNormalOrder();
        }
        
        /// <summary>
        /// Construct LadderOperators from sequence of ladder operators.
        /// </summary>
        /// <param name="setSequence">Sequence of ladder operators.</param>
        public NormalOrderedLadderOperators(List<LadderOperator> setSequence, int setCoefficient = 1) : base(setSequence, setCoefficient)
        {
            ExceptionIfNotInNormalOrder();
        }


        /// <summary>
        /// Construct LadderOperators from sequence of ladder operators.
        /// </summary>
        /// <param name="setSequence">Sequence of ladder operators.</param>
        public NormalOrderedLadderOperators(List<(LadderOperator.Type, int)> set) : base(set)
        {
            ExceptionIfNotInNormalOrder();
        }

        /// <summary>
        /// FermionTerm constructor that assumes normal-ordered fermionic 
        /// creation and annihilation operators, and that the number of
        /// creation an annihilation operators are equal.
        /// </summary>
        public NormalOrderedLadderOperators(IEnumerable<int> indices) : base(indices)
        {
            ExceptionIfNotInNormalOrder();
        }
        #endregion


        /// <summary>
        /// This throws an ArgumentException if the operators in NormalOrderedLadderOperators are not normal-ordered.
        /// </summary>
        private void ExceptionIfNotInNormalOrder()
        {
            if (!base.IsInNormalOrder())
            {
                throw new ArgumentException("NormalOrderedLadderOperators must contatin normal-ordered LadderOperators");
            }
        }
        /// <summary>
        ///  Checks if raising operators indices are in ascending order, 
        ///  then if lowering operator indices are in descending order.
        /// </summary>
        /// <returns>
        /// Returns <c>true</c> this condition is satisfied.
        /// Returns <c>false</c> otherwise.
        /// </returns>
        public bool IsInIndexOrder()
        {
            return IsInIndexCreationCanonicalOrder() && IsInIndexAnnihilationCanonicalOrder();
        }

        /// <summary>
        ///  Checks whether the creation operator sequence of a <see cref="LadderOperators"/> is in 
        ///  canonical order. This means
        ///  1) <c>SpinOrbital</c> is sorted in ascending order for the creation operators.
        /// </summary>
        /// <returns><c>true</c> if the creation opeartor sequence of a <c>LadderOperators</c> is in 
        /// canonical order. <c>false</c> otherwise</returns>
        private bool IsInIndexCreationCanonicalOrder()
        {
            return sequence.Where(o => o.type == LadderOperator.Type.u).Select(o => o.index).IsIntArrayAscending();
        }

        /// <summary>
        ///  Checks whether the annihilation operator sequence of a <see cref="LadderOperators"/> is in 
        ///  canonical order. This means
        ///  1) <c>SpinOrbital</c> is sorted in descending order for the annihilation operators.
        /// </summary>
        /// <returns><c>true</c> if the annihilation opeartor sequence of a <c>LadderOperators</c> is in 
        /// canonical order. <c>false</c> otherwise</returns>
        private bool IsInIndexAnnihilationCanonicalOrder()
        {
            return sequence.Where(o => o.type == LadderOperator.Type.d).Select(o => o.index).Reverse().IsIntArrayAscending();
        }


        public LadderOperators CreateIndexOrder()
        {
            return CreateIndexOrder(this);
        }

        public static NormalOrderedLadderOperators CreateIndexOrder(NormalOrderedLadderOperators ladderOperators)
        {
            var tmp = new NormalOrderedLadderOperators(ladderOperators);
            //var left = sequence.Where(o => o.type == LadderOperator.Type.u).OrderBy(o => o.index);
            //var right = sequence.Where(o => o.type == LadderOperator.Type.d).OrderBy(o => o.index).Reverse();
            //var normalOrdered = new LadderOperators(left.Concat(right));
            //return (0, normalOrdered);

            // Check that LadderOperators is normal-ordered.
            /*
            if (!tmp.IsInNormalOrder())
            {
                throw new System.ArgumentException(
                    $"ToCanonicalOrder() assumes input is normal-ordered. This is currently not satisfied."
                    );
            }*/
            // Check that LadderOperators spin-orbital indices are in canonical order.
            if (!tmp.IsInIndexOrder())
            {

                var upArrayIndices = tmp.sequence.Select((op, idx) => new { op, idx }).Where(x => x.op.type == LadderOperator.Type.u).Select(x => x.idx).ToArray();
                var downArrayIndices = tmp.sequence.Select((op, idx) => new { op, idx }).Where(x => x.op.type == LadderOperator.Type.d).Select(x => x.idx).ToArray();

                // Bubble sort spin-orbital indices of creation operator.
                while (!tmp.IsInIndexCreationCanonicalOrder())
                {
                    for (int idx = 0; idx < upArrayIndices.Count() - 1; idx++)
                    {
                        if (tmp.sequence.ElementAt(upArrayIndices.ElementAt(idx)).index > tmp.sequence.ElementAt(upArrayIndices.ElementAt(idx + 1)).index)
                        {
                            var tmpLadderOperator = tmp.sequence.ElementAt(upArrayIndices.ElementAt(idx));
                            tmp.sequence[upArrayIndices.ElementAt(idx)] = tmp.sequence[upArrayIndices.ElementAt(idx + 1)];
                            tmp.sequence[upArrayIndices.ElementAt(idx + 1)] = tmpLadderOperator;
                            tmp.coefficient = -1 * tmp.coefficient;
                        }
                    }
                }

                // Bubble sort spin-orbital indices of annihilation operator.
                while (!tmp.IsInIndexAnnihilationCanonicalOrder())
                {
                    for (int idx = 0; idx < downArrayIndices.Length - 1; idx++)
                    {
                        if (tmp.sequence.ElementAt(downArrayIndices.ElementAt(idx)).index < tmp.sequence.ElementAt(downArrayIndices.ElementAt(idx + 1)).index)
                        {
                            var tmpLadderOperator = tmp.sequence.ElementAt(downArrayIndices.ElementAt(idx));
                            tmp.sequence[downArrayIndices.ElementAt(idx)] = tmp.sequence[downArrayIndices.ElementAt(idx + 1)];
                            tmp.sequence[downArrayIndices.ElementAt(idx + 1)] = tmpLadderOperator;
                            tmp.coefficient = -1 * tmp.coefficient;
                        }
                    }
                }
            }
            return tmp;
        }
    }

}



