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
    public partial class HermitianFermionTerm : FermionTerm
    {
        /// <summary>
        /// FermionTerm constructor that assumes normal-ordered fermionic 
        /// creation and annihilation operators, and that the number of
        /// creation an annihilation operators are equal.
        /// </summary>
        public HermitianFermionTerm(IEnumerable<int> indices, bool sort = true)
        {
            var length = indices.Count();
            if (length % 2 == 1)
            {
                throw new System.ArgumentException(
                    $"SpinOrbital array of length {length} must be of even length."
                    );
            }
            var tmp = new FermionTerm(indices.Select((o, idx) => new LadderOperator((idx < length / 2 ? LadderOperator.Type.u : LadderOperator.Type.d, o))).ToList());

            sequence = tmp.sequence;
        }

        
        public (int, HermitianFermionTerm) ToCanonicalOrderFromNormalOrder()
        {
            var (coeff, tmp) = base.ToCanonicalOrderFromNormalOrder();

            // Take Hermitian conjugate if still not in canonical order. 
            if (!tmp.IsInCanonicalOrder())
            {
                tmp.sequence = tmp.sequence.Select(o => (o.type == LadderOperator.Type.d ? LadderOperator.Type.u : LadderOperator.Type.d, o.index)).Select(o => new LadderOperator(o)).Reverse().ToList();
            }
            return (coeff, new HermitianFermionTerm(tmp.sequence));
        }

        /// <summary>
        ///  Converts a <c>FermionTerm</c> to canonical order. This generates
        ///  new terms and modifies the coefficient as needed.
        /// </summary>
        public List<(int, HermitianFermionTerm)> ToCanonicalOrder()
        {
            var newTerms = new HermitianFermionTerm(sequence).ToCanonicalOrder();
            // Anti-commutes spin-orbital indices to canonical order
            // and changes the sign of the coefficient as necessay.
            for (int idx = 0; idx < newTerms.Count(); idx++)
            {
                var (coeff, tmp2) = newTerms[idx].Item2.ToCanonicalOrderFromNormalOrder();
                newTerms[idx] = (newTerms[idx].Item1 * coeff, tmp2);
            }
            return newTerms;
        }

    }
}



