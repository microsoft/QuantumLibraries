// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using System;
using System.Linq;
using System.Collections.Generic;

namespace Microsoft.Quantum.Chemistry
{


    public static partial class Extensions
    {
        #region Reordering methods

        /// <summary>
        ///  Converts a <see cref="LadderSequence"/> to normal order. 
        ///  In general, this can generate new terms and modifies the coefficient.
        /// </summary>
        public static HashSet<NormalOrderedLadderSequence> CreateNormalOrder(this LadderSequence ladderOperator)
        {
            // Recursively anti-commute creation to the left.
            var TmpTerms = new Stack<LadderSequence>();
            var NewTerms = new HashSet<NormalOrderedLadderSequence>();

            TmpTerms.Push(new LadderSequence(ladderOperator));

            // Anti-commutes creation and annihilation operators to canonical order
            // and creates new terms if spin-orbital indices match.
            while (TmpTerms.Any())
            {
                var tmpTerm = TmpTerms.Pop();
                if (tmpTerm.IsInNormalOrder())
                {
                    NewTerms.Add(new NormalOrderedLadderSequence(tmpTerm));
                }
                else
                {
                    // Anticommute creation and annihilation operators.
                    for (int i = 0; i < tmpTerm.sequence.Count() - 1; i++)
                    {
                        if ((int)tmpTerm.sequence.ElementAt(i).type > (int)tmpTerm.sequence.ElementAt(i + 1).type)
                        {
                            // If the two elements have the same spin orbital index, generate a new term.
                            if (tmpTerm.sequence.ElementAt(i).index == tmpTerm.sequence.ElementAt(i + 1).index)
                            {
                                var newTerm = new LadderSequence(tmpTerm);
                                newTerm.sequence.RemoveRange(i, 2);
                                TmpTerms.Push(newTerm);
                            }

                            // Swap the two elements and flip sign of the coefficient.
                            var tmpOp = tmpTerm.sequence.ElementAt(i + 1);
                            tmpTerm.sequence[i + 1] = tmpTerm.sequence.ElementAt(i);
                            tmpTerm.sequence[i] = tmpOp;
                            tmpTerm.coefficient *= -1;
                        }
                    }
                    TmpTerms.Push(tmpTerm);
                }
            }
            return NewTerms;
        }

        /// <summary>
        ///  Converts a <see cref="NormalOrderedLadderSequence"/> to index order. 
        ///  In general, this can generate new terms and modifies the coefficient.
        /// </summary>
        public static IndexOrderedLadderSequence CreateIndexOrder(this NormalOrderedLadderSequence ladderOperators)
        {
            var tmp = new IndexOrderedLadderSequence();
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
        ///  Converts a <see cref="LadderSequence"/> to normal order, then index order. 
        ///  In general, this can generate new terms and modifies the coefficient.
        /// </summary>
        public static HashSet<IndexOrderedLadderSequence> CreateIndexOrder(LadderSequence ladderOperator)
        {
            return CreateNormalOrder(ladderOperator).Select(o => CreateIndexOrder(o));
        }
        #endregion
    }

}



