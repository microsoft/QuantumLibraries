// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using System;
using System.Linq;
using System.Collections.Generic;

using Microsoft.Quantum.Chemistry.LadderOperators;
using Microsoft.Quantum.Chemistry;
using Microsoft.Quantum.Chemistry;

namespace Microsoft.Quantum.Chemistry.Fermion
{
    using static Microsoft.Quantum.Chemistry.Extensions;
    using FermionOperator = LadderOperator<int>;

    /// <summary>
    /// Class representing a sequence of fermionic raising and lowering operators, subject to the additional constraints: 
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
    /// <item>
    /// <description>
    /// Hermitian, and is assumed to be implicitly summed with its Hermitian conjugate if not explicitly Hermitian.
    /// </description>
    /// </item>
    /// </summary>
    public class HermitianFermionTerm : FermionTerm, ITermIndex<TermType.Fermion, HermitianFermionTerm>
    {
        #region Constructors
        /// <summary>
        /// Constructor for empty instance.
        /// </summary>
        public HermitianFermionTerm() : base() { }

        /// <summary>
        /// Constructs a Hermitian fermion term from a normal-ordered sequence of ladder operators.
        /// </summary>
        /// <param name="ladderOperators">Normal-ordered sequence of ladder operators.</param>
        public HermitianFermionTerm(LadderSequence<int> ladderOperators) : base(ladderOperators) {
            NormalizeToCanonicalOrder();
        }

        /// <summary>
        /// Constructs a Hermitian fermion term from a normal-ordered sequence of ladder operators.
        /// </summary>
        /// <param name="ladderOperators">Hermitian normal-ordered sequence of ladder operators.</param>
        public HermitianFermionTerm(IEnumerable<FermionOperator> ladderOperators, int setSign = 1) : base(ladderOperators, setSign) {
            NormalizeToCanonicalOrder();
        }

        /// <summary>
        /// Constructs a Hermitian fermion term with an equal number of creation and annihilation operators
        /// from a sequence of integers.
        /// </summary>
        /// <param name="ladderOperators">Hermitian normal-ordered sequence of ladder operators.</param>
        public HermitianFermionTerm(IEnumerable<int> indices) : base(indices.ToLadderSequence())
        {
            NormalizeToCanonicalOrder();
        }

        public static implicit operator HermitianFermionTerm((RaisingLowering, int)[] setSequence)
            => new HermitianFermionTerm(setSequence.Select(o => new FermionOperator(o)));
        #endregion

        /// <summary>
        /// Returns the sign of this fermion term.
        /// </summary>
        /// <returns>The sign of the fermion term.</returns>
        public int Sign => Coefficient;

        /// <summary>
        /// Sets the sign of this fermion term to one.
        /// </summary>
        public void ResetSign()
        {
            Coefficient = 1;
        }

        /// <summary>
        ///  Checks if raising operators indices are in ascending order, 
        ///  then if lowering operator indices are in descending order.
        /// </summary>
        /// <returns>
        /// Returns <c>true</c> this condition is satisfied.
        /// Returns <c>false</c> otherwise.
        /// </returns>
        /// <remarks>
        /// This should always return <c>true</c> when invoked outside this class.
        /// </remarks>
        public bool IsInCanonicalOrder()
        {
            if (IsInIndexOrder())
            {
                var creationSequence = Sequence.Where(o => o.Type == RaisingLowering.u).Select(o => o.Index);
                var annihilationSequence = Sequence.Where(o => o.Type == RaisingLowering.d).Select(o => o.Index);
                if (creationSequence.Count() == annihilationSequence.Count())
                {
                    if (CompareArray(creationSequence, annihilationSequence.Reverse()) > 0)
                    {
                        return false;
                    }
                }
                else if(creationSequence.Count() < annihilationSequence.Count())
                {
                    return false;
                }
                return true;
            }
            return false;
        }

        /// <summary>
        /// Additional sort using the Hermitian conjugate.
        /// </summary>
        private void NormalizeToCanonicalOrder()
        {
            // Take Hermitian Conjugate    
            if (!IsInCanonicalOrder())
            {
                Sequence = Sequence.Select(o => (o.Type == RaisingLowering.d ? RaisingLowering.u : RaisingLowering.d, o.Index)).Select(o => new FermionOperator(o)).Reverse().ToList();
            }
        }

        /// <summary>
        /// Creates a copy of this instance.
        /// </summary>
        public HermitianFermionTerm Clone()
        {
            var newIndices = ToIndices().ToArray().Clone<int>();
            var newOperators = ToRaisingLowering().ToArray().Clone<RaisingLowering>();
            return new HermitianFermionTerm()
            {
                Coefficient = Coefficient,
                Sequence = newOperators.Zip(newIndices, (a, b) => new LadderOperator<int>(a, b)).ToList()
            };
        }

        /// <summary>
        /// Return the category of this fermion term.
        /// </summary>
        /// <returns>Category of fermion term.</returns>
        public TermType.Fermion TermType
        {
            get
            {
                var length = Sequence.Count();
                var uniqueIndices = this.UniqueIndices();

                switch (length)
                {
                    case 0:
                        return Chemistry.TermType.Fermion.Identity;
                    case 2:
                        switch (uniqueIndices)
                        {
                            case 1:
                                return Chemistry.TermType.Fermion.PP;
                            case 2:
                                return Chemistry.TermType.Fermion.PQ;
                            default:
                                throw new ArgumentException("Attempted to classify unknown fermion term.");
                        }
                    case 4:
                        switch (uniqueIndices)
                        {
                            case 2:
                                return Chemistry.TermType.Fermion.PQQP;
                            case 3:
                                return Chemistry.TermType.Fermion.PQQR;
                            case 4:
                                return Chemistry.TermType.Fermion.PQRS;
                            default:
                                throw new ArgumentException("Attempted to classify unknown fermion term.");
                        }
                    default:
                        throw new ArgumentException("Attempted to classify unknown fermion term.");
                }
            }
        }
        
        #region Equality Testing
        public override bool Equals(object obj) => (obj is HermitianFermionTerm x) ? Equals(x) : false;

        public bool Equals(HermitianFermionTerm x) => base.Equals(x);

        public override int GetHashCode() => base.GetHashCode();
        #endregion

    }


}



