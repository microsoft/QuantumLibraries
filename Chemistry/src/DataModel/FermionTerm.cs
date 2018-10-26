// Copyright (c) Microsoft Corporation. All rights reserved. Licensed under the
// Microsoft Software License Terms for Microsoft Quantum Simulation Library (Preview).
// See LICENSE.md in the project root for license information.


using Microsoft.Quantum.Simulation.Core;

using System;
using System.Linq;
using System.Collections.Generic;
using Microsoft.Quantum.Chemistry;

namespace Microsoft.Quantum.Chemistry
{
    /// <summary>
    /// Represents a term in a Fermion Hamiltonian.
    /// </summary>
    [Serializable]
    public struct FermionTerm
    {
        /// <summary>
        /// <c>Int64[] CreationAnnihilation</c> represents a sequence of creation or annihilation operators.
        /// For example, <c>{1,1,0,1,0}</c> represents h * a^\dag_i a^\dag_j a_k a^\dag_l a_m.
        /// </summary>
        public Int64[] CreationAnnihilationIndices;
        /// <summary>
        /// <c>Int64[] SpinOrbital</c> represents the index of a sequence of creation or annihilation operators.
        /// For example, <c>new SpinOrbital[] {i,j,k,l,m}</c> represents the subscript index of, say, h * a^\dag_i a^\dag_j a_k a^\dag_l a_m.
        /// </summary>
        public SpinOrbital[] SpinOrbitalIndices;
        /// <summary>
        /// <c>Double coeff</c> represents the coefficient of a sequence of creation or annihilation operators.
        /// For example, <c>coeff = h</c> represents h * a^\dag_i a^\dag_j a_k a^\dag_l a_m.
        /// </summary>
        public Double coeff;

        /// <summary>
        /// Returns a human-readable description of a spin-orbital.
        /// </summary>
        public override string ToString() =>
            $"([{String.Join(",", CreationAnnihilationIndices)}], " +
            $"[{String.Join(",", SpinOrbitalIndices)}], " +
            $"{coeff})";

        /// <summary>
        /// FermionTerm constructor. 
        /// </summary>
        public FermionTerm(Int64 nOrbitals, Int64[] caArray, Int64[] fermionIdxArray, Double coeffIn)
        {
            CreationAnnihilationIndices = caArray;
            SpinOrbitalIndices = SpinOrbital.ToSpinOrbitals(nOrbitals, fermionIdxArray);
            coeff = coeffIn;
            if (!IsValid())
            {
                throw new System.ArgumentException(
                    $"Invalid FermionTerm specified. Length of caArray {caArray} and soArray {fermionIdxArray} must be equal.",
                    "caArray"
                    );
            }
        }

        /// <summary>
        /// FermionTerm constructor that assumes normal-ordered fermionic 
        /// creation and annihilation operators, and that the number of
        /// creation an annihilation operators are equal.
        /// </summary>
        public FermionTerm(IEnumerable<SpinOrbital> SpinOrbitals, Double coeffIn)
        {
            var length = SpinOrbitals.Count();
            if (length % 2 == 1)
            {
                throw new System.ArgumentException(
                    $"SpinOrbital array of length {length} must be of even length."
                    );
            }
            // Create normal-ordered sequence of operators with
            // an equal number of creation and annihilation operators
            CreationAnnihilationIndices = new Int64[length];
            for (int i = 0; i < length / 2; i++)
            {
                CreationAnnihilationIndices[i] = 1;
                CreationAnnihilationIndices[length - i - 1] = 0;
            }
            SpinOrbitalIndices = SpinOrbitals.ToArray();
            coeff = coeffIn;

            // Anti-commutes spin-orbital indices to canonical order.
            ToSpinOrbitalCanonicalOrder();
        }

        /// <summary>
        /// FermionTerm constructor.
        /// </summary>
        public FermionTerm(Int64[] caArray, SpinOrbital[] soArray, Double coeffIn)
        {
            CreationAnnihilationIndices = caArray;
            SpinOrbitalIndices = soArray;
            coeff = coeffIn;
            if (!IsValid())
            {
                throw new System.ArgumentException(
                    $"Invalid FermionTerm specified. Length of caArray {caArray} and soArray {soArray} must be equal.",
                    "caArray"
                    );
            }
        }

        /// <summary>
        /// Concatenates two Fermion terms.
        /// </summary>
        /// <param name="left">Left Fermion term `x`</param>
        /// <param name="right">Right Fermion term  `y`</param>
        /// <returns>Returns new <see cref="FermionTerm"/> `xy` where coefficients and 
        /// Fermion operators are multipled together.</returns>
        public FermionTerm Concatenate(FermionTerm left, FermionTerm right)
        {
            return new FermionTerm(
                left.CreationAnnihilationIndices.Concat(right.CreationAnnihilationIndices).ToArray(),
                left.SpinOrbitalIndices.Concat(right.SpinOrbitalIndices).ToArray(),
                left.coeff * right.coeff
                );
        }

        /// <summary>
        /// Method for incrementing the coefficient of a <c>FermionTerm</c>.
        /// </summary>
        public void AddCoeff(Double addCoeff)
        {
            coeff += addCoeff;
        }

        public bool IsValid()
        {
            if (CreationAnnihilationIndices.Length == SpinOrbitalIndices.Length)
            {
                return true;
            }
            else
            {
                return false;
            }
        }

        /// <summary>
        /// Count the number of unique <see cref="SpinOrbital"/> indices.
        /// </summary>
        /// <returns>Number of number of unique <see cref="SpinOrbital"/> indices.</returns>
        public Int64 GetUniqueIndices()
        {
            return SpinOrbitalIndices.ToInts().Distinct().Count();
        }

        /// <summary>
        ///  Checks whether a <c>FermionTerm</c> is in canonical order. This means
        ///  1) <c>CreationAnnihilation</c> is a sequence of ones followed by a sequence of zeroes.
        ///  2) <c>SpinOrbital</c> is sorted in ascending order for the creation operators.
        ///  3) <c>SpinOrbital</c> is sorted in descending order for the annihilation operators.
        ///  4) The spin-orbital index of the first creation operator is smaller than
        ///  the spin-orbital index of the last annihilation operator.
        /// </summary>
        /// <returns><c>true</c> if <c>FermionTerm</c> is in canonical order. <c>false</c> otherwise</returns>
        public bool IsInCanonicalOrder()
        {
            if (!IsInNormalOrder())
            {
                return false;
            }
            if (!IsInSpinOrbitalCanonicalOrder())
            {
                return false;
            }
            return true;
        }

        // <summary>
        ///  Checks whether the <see cref="CreationAnnihilationIndices"/> sequence of a <c>FermionTerm</c> is 
        ///  in canonical order. This means
        ///  1) <see cref="CreationAnnihilationIndices"/> is a sequence of ones followed by a sequence of zeroes.
        /// <returns><c>true</c> if the <see cref="CreationAnnihilationIndices"/> sequence of a <c>FermionTerm</c> 
        /// is in canonical order. <c>false</c> otherwise</returns>
        public bool IsInNormalOrder()
        {
            if (CreationAnnihilationIndices.Length == 0)
            {
                return true;
            }
            else if (!CreationAnnihilationIndices.Reverse().IsIntArrayAscending())
            {
                return false;
            }
            return true;
        }

        /// <summary>
        ///  Checks whether the <see cref="SpinOrbitalIndices"/> sequence of a <see cref="FermionTerm"/> is in 
        ///  canonical order. This means
        ///  1) <c>SpinOrbital</c> is sorted in ascending order for the creation operators.
        ///  2) <c>SpinOrbital</c> is sorted in descending order for the annihilation operators.
        ///  3) The spin-orbital index of the first creation operator is smaller than
        ///  the spin-orbital index of the last annihilation operator. If there is a tie,
        ///  look at the next operators.
        /// </summary>
        /// <returns><c>true</c> if the <see cref="SpinOrbitalIndices"/> sequence of a <c>FermionTerm</c> is in 
        /// canonical order. <c>false</c> otherwise</returns>
        public bool IsInSpinOrbitalCanonicalOrder()
        {
            if (!IsInSpinOrbitalCreationCanonicalOrder())
            {
                return false;
            }
            if (!IsInSpinOrbitalAnnihilationCanonicalOrder())
            {
                return false;
            }

            var tmp = CreationAnnihilationIndices.Zip(SpinOrbitalIndices.ToInts(), (a, b) => (a, b));
            var creation = tmp.Where(o => o.a == 1).Select(o => o.b);
            var annihilation = tmp.Where(o => o.a == 0).Select(o => o.b);
            if (creation.Count() == annihilation.Count())
            {
                if (Extensions.CompareIntArray(creation, annihilation.Reverse()) > 0)
                {
                    return false;
                }
            }
            return true;
        }

        /// <summary>
        ///  Checks whether the creation operator sequence of a <see cref="FermionTerm"/> is in 
        ///  canonical order. This means
        ///  1) <c>SpinOrbital</c> is sorted in ascending order for the creation operators.
        /// </summary>
        /// <returns><c>true</c> if the creation opeartor sequence of a <c>FermionTerm</c> is in 
        /// canonical order. <c>false</c> otherwise</returns>
        private bool IsInSpinOrbitalCreationCanonicalOrder()
        {
            var tmp = CreationAnnihilationIndices.Zip(SpinOrbitalIndices.ToInts(), (a, b) => (a, b));
            if (!tmp.Where(o => o.a == 1).Select(o => o.b).IsIntArrayAscending())
            {
                return false;
            }
            return true;
        }

        /// <summary>
        ///  Checks whether the annihilation operator sequence of a <see cref="FermionTerm"/> is in 
        ///  canonical order. This means
        ///  1) <c>SpinOrbital</c> is sorted in descending order for the annihilation operators.
        /// </summary>
        /// <returns><c>true</c> if the annihilation opeartor sequence of a <c>FermionTerm</c> is in 
        /// canonical order. <c>false</c> otherwise</returns>
        private bool IsInSpinOrbitalAnnihilationCanonicalOrder()
        {
            var tmp = CreationAnnihilationIndices.Zip(SpinOrbitalIndices.ToInts(), (a, b) => (a, b));
            if (!tmp.Where(o => o.a == 0).Select(o => o.b).Reverse().IsIntArrayAscending())
            {
                return false;
            }
            return true;
        }

        /// <summary>
        /// Given spin-orbital indices to a normal-ordered fermionic operator, 
        /// this anti-commutes the creation operators and the annihilation orders
        /// so that indices are in canonical-order.
        /// </summary>
        /// <param name="SpinOrbitals">Sequence of spin-orbital indices</param>
        public void ToSpinOrbitalCanonicalOrder()
        {
            // Check that FermionTerm is normal-ordered.
            if (!IsInNormalOrder())
            {
                throw new System.ArgumentException(
                    $"ToSpinOrbitalCanonicalOrder() assumes that Fermionic operators in FermionTerm are normal-ordered." +
                    $"This is currently not satisfied."
                    );
            }

            // Check that FermionTerm spin-orbital indices are in canonical order.
            if (!IsInSpinOrbitalCanonicalOrder())
            {

                var CreationIndex = CreationAnnihilationIndices.Select((op, idx) => new { op, idx }).Where(x => x.op == 1).Select(x => x.idx).ToArray();
                var AnnihilationIndex = CreationAnnihilationIndices.Select((op, idx) => new { op, idx }).Where(x => x.op == 0).Select(x => x.idx).ToArray();

                // Bubble sort spin-orbital indices of creation operator.
                while (!IsInSpinOrbitalCreationCanonicalOrder())
                {
                    for (int idx = 0; idx < CreationIndex.Length - 1; idx++)
                    {
                        if (SpinOrbitalIndices.ElementAt(CreationIndex.ElementAt(idx)).ToInt() > SpinOrbitalIndices.ElementAt(CreationIndex.ElementAt(idx + 1)).ToInt())
                        {
                            var tmpSpinOrbital = SpinOrbitalIndices.ElementAt(CreationIndex.ElementAt(idx));
                            SpinOrbitalIndices[CreationIndex.ElementAt(idx)] = SpinOrbitalIndices[CreationIndex.ElementAt(idx + 1)];
                            SpinOrbitalIndices[CreationIndex.ElementAt(idx + 1)] = tmpSpinOrbital;
                            coeff = -1.0 * coeff;
                        }
                    }
                }

                // Bubble sort spin-orbital indices of annihilation operator.
                while (!IsInSpinOrbitalAnnihilationCanonicalOrder())
                {
                    for (int idx = 0; idx < AnnihilationIndex.Length - 1; idx++)
                    {
                        if (SpinOrbitalIndices.ElementAt(AnnihilationIndex.ElementAt(idx)).ToInt() < SpinOrbitalIndices.ElementAt(AnnihilationIndex.ElementAt(idx + 1)).ToInt())
                        {
                            var tmpSpinOrbital = SpinOrbitalIndices.ElementAt(AnnihilationIndex.ElementAt(idx));
                            SpinOrbitalIndices[AnnihilationIndex.ElementAt(idx)] = SpinOrbitalIndices[AnnihilationIndex.ElementAt(idx + 1)];
                            SpinOrbitalIndices[AnnihilationIndex.ElementAt(idx + 1)] = tmpSpinOrbital;
                            coeff = -1.0 * coeff;
                        }
                    }
                }

                // Take Hermitian conjugate if still not in canonical order. 
                if (!IsInSpinOrbitalCanonicalOrder()) { 
                    CreationAnnihilationIndices = CreationAnnihilationIndices.Reverse().Select(o => (o + 1) % 2).ToArray();
                    SpinOrbitalIndices = SpinOrbitalIndices.Reverse().ToArray();
                }
            }
        }

        /// <summary>
        ///  Converts a <c>FermionTerm</c> to canonical order. This generates
        ///  new terms and modifies the coefficient as needed.
        /// </summary>
        public List<FermionTerm> ToCanonicalOrder()
        {
            // Step 1: anti-commute creation to the left.
            // Step 2: sort to canonical order

            var TmpTerms = new Stack<FermionTerm>();
            var NewTerms = new List<FermionTerm>();

            TmpTerms.Push(this);

            // Anti-commutes creation and annihilation operators to canonical order
            // and creates new terms if spin-orbital indices match.
            while (TmpTerms.Any())
            {
                var tmpTerm = TmpTerms.Pop();
                if (tmpTerm.IsInNormalOrder())
                {
                    NewTerms.Add(tmpTerm);
                }
                else
                {
                    // Anticommute creation and annihilation operators.
                    for (int i = 0; i < tmpTerm.CreationAnnihilationIndices.Count() - 1; i++)
                    {
                        if (tmpTerm.CreationAnnihilationIndices.ElementAt(i) < tmpTerm.CreationAnnihilationIndices.ElementAt(i + 1))
                        {
                            var antiCommutedCreationAnnihilationIndices = tmpTerm.CreationAnnihilationIndices.ToList();
                            var antiCommutedSpinOrbitalIndices = tmpTerm.SpinOrbitalIndices.ToList();
                            // Swap the two elements and flip sign of the coefficient.
                            antiCommutedCreationAnnihilationIndices[i + 1] = tmpTerm.CreationAnnihilationIndices.ElementAt(i);
                            antiCommutedCreationAnnihilationIndices[i] = tmpTerm.CreationAnnihilationIndices.ElementAt(i + 1);
                            antiCommutedSpinOrbitalIndices[i + 1] = tmpTerm.SpinOrbitalIndices.ElementAt(i);
                            antiCommutedSpinOrbitalIndices[i] = tmpTerm.SpinOrbitalIndices.ElementAt(i + 1);

                            var antiCommutedTerm = new FermionTerm(antiCommutedCreationAnnihilationIndices.ToArray(), antiCommutedSpinOrbitalIndices.ToArray(), -1.0 * tmpTerm.coeff);

                            TmpTerms.Push(antiCommutedTerm);

                            // If the two elements have the same spin orbital index, generate a new term.
                            if (antiCommutedSpinOrbitalIndices.ElementAt(i).Equals(antiCommutedSpinOrbitalIndices.ElementAt(i + 1)))
                            {
                                var newCreationAnnihilationIndices = antiCommutedCreationAnnihilationIndices.ToList();
                                var newSpinOrbitalIndices = antiCommutedSpinOrbitalIndices.ToList();
                                newCreationAnnihilationIndices.RemoveRange(i, 2);
                                newSpinOrbitalIndices.RemoveRange(i, 2);

                                var newTerm = new FermionTerm(newCreationAnnihilationIndices.ToArray(), newSpinOrbitalIndices.ToArray(), tmpTerm.coeff);

                                TmpTerms.Push(newTerm);
                            }
                            break;
                        }
                    }
                }
            }

            // Anti-commutes spin-orbital indices to canonical order
            // and changes the sign of the coefficient as necessay.
            for (int idx = 0; idx < NewTerms.Count(); idx++)
            {
                var tmp = NewTerms[idx];
                tmp.ToSpinOrbitalCanonicalOrder();
                NewTerms[idx] = tmp;
            }
            return NewTerms;
        }

        /// <summary>
        /// Gets the <see cref="FermionTermType"/> of this term.
        /// </summary>
        /// <returns>The <see cref="FermionTermType"/> of this term.</returns>
        public FermionTermType GetFermionTermType()
        {
            // Convets term to canonical order. 
            ToSpinOrbitalCanonicalOrder();

            //
            return new FermionTermType(GetUniqueIndices(), CreationAnnihilationIndices);
        }

        /// <summary>
        /// Boolean equality operator definition.
        /// </summary>
        public static bool operator ==(FermionTerm x, FermionTerm y)
        {
            return x.coeff == y.coeff
                && x.CreationAnnihilationIndices.SequenceEqual(y.CreationAnnihilationIndices)
                && x.SpinOrbitalIndices.SequenceEqual(y.SpinOrbitalIndices);
        }
        /// <summary>
        /// Boolean inequality operator definition.
        /// </summary>
        public static bool operator !=(FermionTerm x, FermionTerm y)
        {
            return !(x == y);
        }

        public override bool Equals(object x)
        {
            return Equals((FermionTerm)x);
        }

        public bool Equals(FermionTerm x)
        {
            if (ReferenceEquals(null, x))
            {
                return false;
            }
            else if (ReferenceEquals(this, x))
            {
                return true;
            }
            else if (GetType() != x.GetType())
            {
                return false;
            }
            else
                return this == x;
        }

        public override int GetHashCode()
        {
            return coeff.GetHashCode() ^ CreationAnnihilationIndices.GetHashCode() ^ SpinOrbitalIndices.GetHashCode();
        }
    }

    /// <summary>
    /// Represents the type of a term in a Fermion Hamiltonian.
    /// </summary>
    [Serializable]
    public struct FermionTermType
    {

        /// <summary>
        /// First parameter <c>Int64</c> is the the number of different spin-orbitals.
        /// Second parameter <c>QArray<Int64></c> represents a sequence of creation and annihilation operators.
        /// Third parameter <c>string</c> annotates the term type.
        /// </summary>
        public readonly (Int64, Int64[]) type;

        /// <summary>
        /// Constructor of <c>FermionTermTyp</c>
        /// </summary>
        /// <param name="uniqueTerms">Number of distint indices.</param>
        /// <param name="creationAnnihilation">Sequence of creation and annihilation operators.</param>
        /// <param name="annotation">Text annotation of term type.</param>
        public FermionTermType(Int64 uniqueTerms, Int64[] creationAnnihilation)
        {
            type = (uniqueTerms, creationAnnihilation);
        }

        /// <summary>
        /// Returns the sequence of creation and annihilation operators.
        /// </summary>
        /// <returns>Sequence of creation and annihilation operators of <see cref="FermionTermType></returns>
        public Int64[] GetConjugateSequence()
        {
            return type.Item2;
        }

        /// <summary>
        /// Commonly used term types in a Fermion Hamiltonian.
        /// </summary>
        public sealed class Common
        {
            #region Common term types
            public readonly static FermionTermType IdentityTermType = new FermionTermType(0, new Int64[] { });
            public readonly static FermionTermType PPTermType = new FermionTermType(1, new Int64[] { 1, 0 });
            public readonly static FermionTermType PQTermType = new FermionTermType(2, new Int64[] { 1, 0 });
            public readonly static FermionTermType PQQPTermType = new FermionTermType(2, new Int64[] { 1, 1, 0, 0 });
            public readonly static FermionTermType PQQRTermType = new FermionTermType(3, new Int64[] { 1, 1, 0, 0 });
            public readonly static FermionTermType PQRSTermType = new FermionTermType(4, new Int64[] { 1, 1, 0, 0 });
            #endregion

        }

        /// <summary>
        /// Checks whether a <c>FermionTermType</c> is in canonical order. This means
        /// that <c>CreationAnnihilation</c> is a sequence of all ones followed by a sequence of all zeroes.
        /// </summary>
        /// <returns>Returns <c>true</c> if in canonical oder, and <c>false</c> otherwise.</returns>
        public bool IsInCanonicalOrder()
        {
            if (type.Item1 > type.Item2.Length)
            {
                return false;
            }
            for (int i = 0; i < type.Item2.Length - 1; i++)
            {
                if (type.Item2[i] < type.Item2[i + 1])
                {
                    return false;
                }
            }
            return true;
        }


        /// <summary>
        /// String representation of <see cref="FermionTermType"/>.
        /// </summary>
        /// <returns>String representation of <see cref="FermionTermType"/>.</returns>
        public override string ToString()
        {
            return $"({type.Item1}, [{String.Join(",", type.Item2)}])";
        }

        /// <summary>
        /// Boolean tests.
        /// </summary>
        public static bool operator ==(FermionTermType x, FermionTermType y)
        {
            return x.type.Item1 == y.type.Item1 && x.type.Item2.SequenceEqual(y.type.Item2);
        }
        public static bool operator !=(FermionTermType x, FermionTermType y)
        {
            return !(x == y);
        }
        public override bool Equals(object x)
        {
            return Equals((FermionTermType)x);
        }

        public bool Equals(FermionTermType x)
        {
            if (ReferenceEquals(null, x))
            {
                return false;
            }
            else if (ReferenceEquals(this, x))
            {
                return true;
            }
            else if (GetType() != x.GetType())
            {
                return false;
            }
            else
                return this == x;
        }

        public override int GetHashCode()
        {
            int h = 19;
            foreach (var i in type.Item2)
            {
                h = h * 53 + (i * type.Item1).GetHashCode();
            }
            return h;
        }
    }
}

namespace Microsoft.Quantum.Chemistry.Comparers
{ 

    /// <summary>
    /// IComparer for <c>FermionTerm</c>.
    /// </summary>
    public class FermionTermIComparer : IComparer<FermionTerm>
    {
        private int nOnes = 0;
        public FermionTermIComparer(int nOnesIn)
        {
            nOnes = nOnesIn;
        }

        public int Compare(FermionTerm x, FermionTerm y)
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
    public class FermionTermComparer : IEqualityComparer<FermionTerm>
    {
        public bool Equals(FermionTerm x, FermionTerm y)
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
        public int GetHashCode(FermionTerm x)
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
    }
}



