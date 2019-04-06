// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using System;
using System.Linq;
using System.Collections.Generic;

namespace Microsoft.Quantum.Chemistry
{
    
    public class IndexOrderedLadderOperators : NormalOrderedLadderOperators
    {
        #region Constructors
        /// <summary>
        /// Constructor for empty ladder operator sequence.
        /// </summary>
        internal IndexOrderedLadderOperators() : base() { }

        /// <summary>
        /// Construct LadderOperators from another LadderOperators.
        /// </summary>
        /// <param name="setSequence">Sequence of ladder operators.</param>
        public IndexOrderedLadderOperators(IndexOrderedLadderOperators ladderOperators)
        {
            // All constructions are pass by value.
            sequence = ladderOperators.sequence.Select(o => o).ToList();
            coefficient = ladderOperators.coefficient;
        }

        /// <summary>
        /// Construct LadderOperators from sequence of ladder operators.
        /// </summary>
        /// <param name="setSequence">Sequence of ladder operators.</param>
        public IndexOrderedLadderOperators(LadderOperatorSequence ladderOperators) : base(ladderOperators)
        {
            ToIndexOrder();
        }
        
        /// <summary>
        /// Construct LadderOperators from sequence of ladder operators.
        /// </summary>
        /// <param name="setSequence">Sequence of ladder operators.</param>
        public IndexOrderedLadderOperators(IEnumerable<LadderOperator> setSequence, int setCoefficient = 1) : base(setSequence, setCoefficient)
        {
            ToIndexOrder();
        }


        /// <summary>
        /// Construct LadderOperators from sequence of ladder operators.
        /// </summary>
        /// <param name="setSequence">Sequence of ladder operators.</param>
        public IndexOrderedLadderOperators(IEnumerable<(LadderOperator.Type, int)> set) : base(set)
        {
            ToIndexOrder();
        }

        /// <summary>
        /// FermionTerm constructor that assumes normal-ordered fermionic 
        /// creation and annihilation operators, and that the number of
        /// creation an annihilation operators are equal.
        /// </summary>
        public IndexOrderedLadderOperators(IEnumerable<int> indices) : base(indices)
        {
            ToIndexOrder();
        }
        #endregion




        #region Reordering methods
        private void ToIndexOrder()
        {
            var ladderTerm = CreateIndexOrder();
            sequence = ladderTerm.sequence;
            coefficient = ladderTerm.coefficient;        
        }

        public IndexOrderedLadderOperators CreateIndexOrder()
        {
            return CreateIndexOrder(this);
        }

        public static IndexOrderedLadderOperators CreateIndexOrder(NormalOrderedLadderOperators ladderOperators)
        {
            var tmp = new IndexOrderedLadderOperators();
            tmp.sequence = ladderOperators.sequence.Select(o => o).ToList();
            tmp.coefficient = ladderOperators.coefficient;
            if (!tmp.IsInIndexOrder())
            {
                var left = tmp.sequence.Select((op, idx) => new { op, idx }).Where(x => x.op.type == LadderOperator.Type.u);
                var right = tmp.sequence.Select((op, idx) => new { op, idx }).Where(x => x.op.type == LadderOperator.Type.d);

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

        /// <summary>
        ///  Converts a <see cref="LadderOperatorSequence"/> to normal order, then index order. 
        ///  In general, this can generate new terms and modifies the coefficient.
        /// </summary>
        public static IEnumerable<IndexOrderedLadderOperators> CreateIndexOrder(LadderOperatorSequence ladderOperator)
        {
            return CreateNormalOrder(ladderOperator).Select(o => CreateIndexOrder(o));
        }
        #endregion
        
    }

}



