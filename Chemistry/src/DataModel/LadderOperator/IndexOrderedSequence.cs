// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.


using System;
using System.Linq;
using System.Collections.Generic;

namespace Microsoft.Quantum.Chemistry.LadderOperators
{

    /// <summary>
    /// Class representing a sequence of raising and lowering operators, subject to the additional constraints: 
    /// <list type="number">
    /// <item>
    /// <description>
    /// Normal-ordered, where all raising operators are to the left of all lowering operators.
    /// </description>
    /// </item>
    /// <item>
    /// <description>
    /// Index-ordered, where are raising(lowering) operators are in ascending(descending) order.
    /// </description>
    /// </item>
    /// </summary>
    public class IndexOrderedSequence<TIndex> : NormalOrderedSequence<TIndex>
        where TIndex : IEquatable<TIndex>, IComparable<TIndex>
    {
        #region Constructors
        /// <summary>
        /// Constructor for empty ladder operator sequence.
        /// </summary>
        public IndexOrderedSequence() : base() { }

        /// <summary>
        /// Constructs an instance from a normal-ordered sequence of ladder operators.
        /// </summary>
        /// <param name="ladderOperators">Normal-ordered sequence of ladder operators.</param>
        public IndexOrderedSequence(LadderSequence<TIndex> ladderOperators) : base(ladderOperators)
        {
            NormalizeToIndexOrder();
        }

        /// <summary>
        /// Constructs an instance from a normal-ordered sequence of ladder operators.
        /// </summary>
        /// <param name="ladderOperators">Normal-ordered sequence of ladder operators.</param>
        public IndexOrderedSequence(IEnumerable<LadderOperator<TIndex>> ladderOperators, int setSign = 1) : base(ladderOperators, setSign)
        {
            NormalizeToIndexOrder();
        }
        #endregion

        #region Reordering methods
        /// <summary>
        ///  Converts a <see cref="NormalOrderedLadderSequence"/> to index order. 
        ///  In general, this can generate new terms and modifies the coefficient.
        /// </summary>
        public void NormalizeToIndexOrder()
        {
            var tmp = new NormalOrderedSequence<TIndex>(this);
            if (!tmp.IsInIndexOrder())
            {
                var left = tmp.Sequence.Select((op, idx) => new { op, idx }).Where(x => x.op.Type == RaisingLowering.u);
                var right = tmp.Sequence.Select((op, idx) => new { op, idx }).Where(x => x.op.Type == RaisingLowering.d);

                var upArrayIndices = tmp.Sequence.Select((op, idx) => new { op, idx }).Where(x => x.op.Type == RaisingLowering.u).Select(x => x.idx).ToArray();
                var downArrayIndices = tmp.Sequence.Select((op, idx) => new { op, idx }).Where(x => x.op.Type == RaisingLowering.d).Select(x => x.idx).ToArray();

                // Bubble sort spin-orbital indices of creation operator.
                while (!tmp.IsInIndexCreationCanonicalOrder())
                {
                    for (int idx = 0; idx < upArrayIndices.Count() - 1; idx++)
                    {
                        if (tmp.Sequence.ElementAt(upArrayIndices.ElementAt(idx)).Index.CompareTo(tmp.Sequence.ElementAt(upArrayIndices.ElementAt(idx + 1)).Index) > 0)
                        {
                            var tmpLadderOperator = tmp.Sequence.ElementAt(upArrayIndices.ElementAt(idx));
                            tmp.Sequence[upArrayIndices.ElementAt(idx)] = tmp.Sequence[upArrayIndices.ElementAt(idx + 1)];
                            tmp.Sequence[upArrayIndices.ElementAt(idx + 1)] = tmpLadderOperator;
                            tmp.Coefficient = -1 * tmp.Coefficient;
                        }
                    }
                }

                // Bubble sort spin-orbital indices of annihilation operator.
                while (!tmp.IsInIndexAnnihilationCanonicalOrder())
                {
                    for (int idx = 0; idx < downArrayIndices.Length - 1; idx++)
                    {
                        if (tmp.Sequence.ElementAt(downArrayIndices.ElementAt(idx)).Index.CompareTo(tmp.Sequence.ElementAt(downArrayIndices.ElementAt(idx + 1)).Index) < 0)
                        {
                            var tmpLadderOperator = tmp.Sequence.ElementAt(downArrayIndices.ElementAt(idx));
                            tmp.Sequence[downArrayIndices.ElementAt(idx)] = tmp.Sequence[downArrayIndices.ElementAt(idx + 1)];
                            tmp.Sequence[downArrayIndices.ElementAt(idx + 1)] = tmpLadderOperator;
                            tmp.Coefficient = -1 * tmp.Coefficient;
                        }
                    }
                }
            }
            Sequence = tmp.Sequence;
            Coefficient = tmp.Coefficient;
        }
        #endregion

    }

}



