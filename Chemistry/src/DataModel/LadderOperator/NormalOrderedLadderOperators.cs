// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using System;
using System.Linq;
using System.Collections.Generic;

namespace Microsoft.Quantum.Chemistry
{
    
    public class NormalOrderedLadderOperators : LadderOperatorSequence
    {
        #region Constructors
        /// <summary>
        /// Constructor for empty ladder operator sequence.
        /// </summary>
        internal NormalOrderedLadderOperators() : base() { }

        /// <summary>
        /// Construct LadderOperators from another LadderOperators.
        /// </summary>
        /// <param name="setSequence">Sequence of ladder operators.</param>
        public NormalOrderedLadderOperators(NormalOrderedLadderOperators ladderOperators)
        {
            // All constructions are pass by value.
            sequence = ladderOperators.sequence.Select(o => o).ToList();
            coefficient = ladderOperators.coefficient;
        }

        /// <summary>
        /// Construct LadderOperators from sequence of ladder operators.
        /// </summary>
        /// <param name="setSequence">Sequence of ladder operators.</param>
        public NormalOrderedLadderOperators(LadderOperatorSequence ladderOperators) : base(ladderOperators)
        {
            ExceptionIfNotInNormalOrder();
        }
        
        /// <summary>
        /// Construct LadderOperators from sequence of ladder operators.
        /// </summary>
        /// <param name="setSequence">Sequence of ladder operators.</param>
        public NormalOrderedLadderOperators(IEnumerable<LadderOperator> setSequence, int setCoefficient = 1) : base(setSequence, setCoefficient)
        {
            ExceptionIfNotInNormalOrder();
        }


        /// <summary>
        /// Construct LadderOperators from sequence of ladder operators.
        /// </summary>
        /// <param name="set">Sequence of ladder operators.</param>
        public NormalOrderedLadderOperators(IEnumerable<(LadderOperator.Type, int)> set) : base(set)
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

        #region Ordering testers
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
        ///  Checks whether the creation operator sequence of a <see cref="LadderOperatorSequence"/> is in 
        ///  canonical order. This means
        ///  1) <c>SpinOrbital</c> is sorted in ascending order for the creation operators.
        /// </summary>
        /// <returns><c>true</c> if the creation opeartor sequence of a <c>LadderOperators</c> is in 
        /// canonical order. <c>false</c> otherwise</returns>
        public bool IsInIndexCreationCanonicalOrder()
        {
            return sequence.Where(o => o.type == LadderOperator.Type.u).Select(o => o.index).IsIntArrayAscending();
        }

        /// <summary>
        ///  Checks whether the annihilation operator sequence of a <see cref="LadderOperatorSequence"/> is in 
        ///  canonical order. This means
        ///  1) <c>SpinOrbital</c> is sorted in descending order for the annihilation operators.
        /// </summary>
        /// <returns><c>true</c> if the annihilation opeartor sequence of a <c>LadderOperators</c> is in 
        /// canonical order. <c>false</c> otherwise</returns>
        public bool IsInIndexAnnihilationCanonicalOrder()
        {
            return sequence.Where(o => o.type == LadderOperator.Type.d).Select(o => o.index).Reverse().IsIntArrayAscending();
        }
        #endregion

        #region Reordering methods

        /// <summary>
        ///  Converts a <see cref="LadderOperatorSequence"/> to normal order. 
        ///  In general, this can generate new terms and modifies the coefficient.
        /// </summary>
        public static IEnumerable<NormalOrderedLadderOperators> CreateNormalOrder(LadderOperatorSequence ladderOperator)
        {
            // Recursively anti-commute creation to the left.
            var TmpTerms = new Stack<LadderOperatorSequence>();
            var NewTerms = new HashSet<NormalOrderedLadderOperators>();
            
            TmpTerms.Push(new LadderOperatorSequence(ladderOperator));

            // Anti-commutes creation and annihilation operators to canonical order
            // and creates new terms if spin-orbital indices match.
            while (TmpTerms.Any())
            {
                var tmpTerm = TmpTerms.Pop();
                if (tmpTerm.IsInNormalOrder())
                {
                    NewTerms.Add(new NormalOrderedLadderOperators(tmpTerm));
                }
                else
                {
                    // Anticommute creation and annihilation operators.
                    for (int i = 0; i < tmpTerm.sequence.Count() - 1; i++)
                    {
                        if ((int) tmpTerm.sequence.ElementAt(i).type > (int) tmpTerm.sequence.ElementAt(i + 1).type)
                        {
                            // If the two elements have the same spin orbital index, generate a new term.
                            if (tmpTerm.sequence.ElementAt(i).index == tmpTerm.sequence.ElementAt(i + 1).index)
                            {
                                var newTerm = new LadderOperatorSequence(tmpTerm);
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
    }

}



