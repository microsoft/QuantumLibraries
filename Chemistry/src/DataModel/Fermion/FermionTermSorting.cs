// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using Microsoft.Quantum.Simulation.Core;

using System;
using System.Linq;
using System.Collections.Generic;
using Microsoft.Quantum.Chemistry;

using Microsoft.Quantum.Chemistry.LadderOperators;

namespace Microsoft.Quantum.Chemistry.Fermion
{ 
    /*
    /// <summary>
    /// IComparer for <c>FermionTerm</c>.
    /// </summary>
    public class FermionTermIComparer : IComparer<FermionTermDeprecated>
    {
        private int nOnes = 0;
        public FermionTermIComparer(int nOnesIn)
        {
            nOnes = nOnesIn;
        }

        public int Compare(FermionTermDeprecated x, FermionTermDeprecated y)
        {
            return Extensions.CompareIntArray(
                x.SpinOrbitalIndices.Take(nOnes).Concat(x.SpinOrbitalIndices.Skip(nOnes).Reverse()).ToInts(),
                y.SpinOrbitalIndices.Take(nOnes).Concat(y.SpinOrbitalIndices.Skip(nOnes).Reverse()).ToInts()
                );
        }
    }

    /// <summary>
    /// Equality comparer for <c>FermionTermType</c>. Two <c>FermionTermType</c>s
    /// are identical when the number of distinct spin-orbit indices are identical
    /// and the sequence of creation and annihilation operators are identical.
    /// </summary>
    [Serializable]
    public class FermionTermTypeComparer : IEqualityComparer<FermionTermType>
    {
        public bool Equals(FermionTermType x, FermionTermType y)
        {
            if (x.type.Item2.SequenceEqual(y.type.Item2) && x.type.Item1 == y.type.Item1)
            {
                return true;
            }
            else
            {
                return false;
            }

        }
        public int GetHashCode(FermionTermType x)
        {
            int h = 19;
            foreach (var i in x.type.Item2)
            {
                h = h * 53 + (i * x.type.Item1).GetHashCode();
            }
            return h;
        }
    }


    /// <summary>
    /// Equality comparer for <c>FermionTerm</c>. Two <c>FermionTerm</c>s
    /// are identical when 1) the sequence of creation and annihilation operators are identical.
    /// 2) The sequence of spin-orbitals are identical.
    /// 3) The coefficients are identical.
    /// </summary>
    public class FermionTermComparer : IEqualityComparer<FermionTermDeprecated>
    {
        public bool Equals(FermionTermDeprecated x, FermionTermDeprecated y)
        {
            if (x.CreationAnnihilationIndices.SequenceEqual(y.CreationAnnihilationIndices) && x.SpinOrbitalIndices.SequenceEqual(y.SpinOrbitalIndices) && x.coeff == y.coeff)
            {
                return true;
            }
            else
            {
                return false;
            }
        }
        public int GetHashCode(FermionTermDeprecated x)
        {
            int h = 19;
            foreach (var i in x.CreationAnnihilationIndices.Select(o => o.GetHashCode()))
            {
                h = h * 31 + i;
            }
            foreach (var i in x.SpinOrbitalIndices.Select(o => o.orbital.GetHashCode()))
            {
                h = h * 17 + i;
            }
            foreach (var i in x.SpinOrbitalIndices.Select(o => o.spin.GetHashCode()))
            {
                h = h * 53 + i;
            }
            return h;
        }
    }*/
}



