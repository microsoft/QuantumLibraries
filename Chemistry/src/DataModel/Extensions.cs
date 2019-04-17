// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using System;
using System.Collections.Generic;
using System.Linq;
using Microsoft.Quantum.Simulation.Core;


namespace Microsoft.Quantum.Chemistry
{
    using TermArray = List<HTerm>;

    public static partial class Extensions
    {
        /// <summary>
        ///      Given a value of an enumeration type, and an action for each
        ///      possible value of that enumeration type, performs the action
        ///      corresponding to the given value.
        /// </summary>
        public static void Map<E>
            (this E @enum, params (E, Action)[] actions)
            where E : struct, IConvertible
        {
            foreach (var (value, action) in actions)
            {
                if (@enum.Equals(value))
                {
                    action();
                    return;
                }
            }
        }

        /// <summary>
        ///      Given a value of an enumeration type, and an function for each
        ///      possible value of that enumeration type, returns the value
        //       returned by the function corresponding to the given value.
        /// </summary>
        public static T Map<E, T>
            (this E @enum, params (E, Func<T>)[] actions)
            where E: struct, IConvertible
        {
            foreach (var (value, action) in actions)
            {
                if (@enum.Equals(value))
                {
                    return action();
                }
            }

            throw new ArgumentException($"Expected {@enum} to be a member of {@enum.GetType()}.");
        }

        #region Extension methods
        /// <summary>
        /// Converts an array of <c>SpinOrbital</c>s into an array of integers representing each spin orbital.
        /// </summary>
        public static Int64[] ToInts(this IEnumerable<SpinOrbital> spinOrbitals, Int64 nOrbitals)
        {
            return spinOrbitals.Select(x => x.ToInt(nOrbitals)).ToArray();
        }

        /// <summary>
        /// Converts an array of <c>SpinOrbital</c>s into an array of integers representing each spin orbital.
        /// </summary>
        internal static Int64[] ToInts(this IEnumerable<SpinOrbital> spinOrbitals)
        {
            return spinOrbitals.Select(x => x.ToInt()).ToArray();
        }

        /// <summary>
        /// Converts an array of (orbital index, spin index) into an array of spin-orbitals.
        /// </summary>
        public static SpinOrbital[] ToSpinOrbitals(this IEnumerable<(int, Spin)> spinOrbitalIndices)
        {
            return spinOrbitalIndices.Select(o => new SpinOrbital(o)).ToArray();
        }
        /// <summary>
        /// Converts an array of (orbital index, spin index) into an array of spin-orbitals.
        /// </summary>
        public static SpinOrbital[] ToSpinOrbitals(this IEnumerable<(int, int)> spinOrbitalIndices)
        {
            return spinOrbitalIndices.Select(o => new SpinOrbital(o)).ToArray();
        }

        /// <summary>
        /// Enumerates all spin-orbitals described by an array of <see cref="OrbitalIntegral"/> by
        /// applying <see cref="OrbitalIntegral.EnumerateSpinOrbitals"/> to each.
        /// </summary>
        /// <param name="orbitalIntegrals">Array of orbital integrals.</param>
        /// <returns>Array of Array of spin-orbitals.</returns>
        public static SpinOrbital[][] EnumerateSpinOrbitals(this IEnumerable<OrbitalIntegral> orbitalIntegrals)
        {
            return orbitalIntegrals.SelectMany(o => o.EnumerateSpinOrbitals()).ToArray();
        }


        /// <summary>
        /// String representation of elements of an <see cref="IEnumerable{T}"/> of <see cref="Int64"/>.
        /// </summary>
        /// <param name="ints"><see cref="IEnumerable{T}"/> of <see cref="Int64"/>.</param>
        /// <returns>String representation of input elements.</returns>
        public static string Print(this IEnumerable<Int64> ints)
        {
            return "[" + string.Join(", ", ints) + "]";
        }

        /// <summary>
        /// Checks whether a sequence of integers {x,y,z,...} is sorted in
        /// ascending order x <= y <= z <= ...
        /// </summary>
        /// <param name="x">Sequence of integers to be checked.</param>
        /// <returns>
        /// <c>true</c> if <paramref name="x"/> is sorted in ascending
        /// order, and <c>false</c> otherwise.
        /// </returns>
        public static bool IsIntArrayAscending(this IEnumerable<Int64> x)
        {
            // return true for empty list.
            if (!x.Any())
            {
                return true;
            }
            var curr = x.First();
            foreach (var next in x.Skip(1))
            {
                if (next < curr)
                {
                    return false;
                }
                curr = next;
            }
            return true;
        }

        #region Accumulators

        /// <summary>
        /// This combines the coefficients of adjacent identical terms in a 
        /// sequence of type <c>(QArray<Int64>,QArray<Double>)</c>. Adjacent terms are
        /// considered identical if their <c>QArray<Int64></c> sequences are identical.
        /// </summary>
        /// <param name="term">Sequence of terms in a Hamiltonian.</param>
        /// <returns>Sequence of terms in a Hamiltonian with adjacent entries of the 
        /// same type merged.</returns>
        public static TermArray AccumulateTermArray(this IEnumerable<HTerm> term)
        {
            TermArray output = new TermArray();
            if (term.Any())
            {
                output.Add(term.First());
                var nElements = 1;
                foreach (var next in term.Skip(1))
                {
                    HTerm curr = output[nElements - 1];

                    if (Enumerable.SequenceEqual(curr.Item1, next.Item1))
                    {
                        var (pqrsSorted, coeffs) = next;
                        for (var idxCoeff = 0; idxCoeff < coeffs.Length; idxCoeff++)
                        {
                            curr = new HTerm((curr.Item1, new QArray<double>(curr.Item2).SetItem(idxCoeff, curr.Item2[idxCoeff] + coeffs[idxCoeff])));
                            output[nElements - 1] = curr;
                        }
                    }
                    else
                    {
                        output.Add(next);
                        nElements++;
                    }
                }
            }
            return output;
        }



        #endregion
        /// <summary>
        /// This combines the coefficients of adjacent identical terms in a 
        /// sequence of type <c>FermionTerm</c> that are sorted in canonical order.
        /// <see cref="FermionTerm.CanonicalSort"/>. Adjacent terms are
        /// considered identical if their creation-annihilation and
        /// spin-orbit sequences are identical.
        /// </summary>
        /// <param name="term">Sequence of terms in a Hamiltonian.</param>
        /// <returns>
        /// Sequence of terms in a Hamiltonian with adjacent entries of the 
        /// same type merged.
        /// </returns>
        public static List<FermionTerm> AccumulateFermionTerm(this IEnumerable<FermionTerm> term)
        {
            var output = new List<FermionTerm>();
            if (term.Any())
            {
                output.Add(term.First());
                var nElements = 1;
                foreach (var next in term.Skip(1))
                {
                    var curr = output[nElements - 1];

                    if (Enumerable.SequenceEqual(curr.SpinOrbitalIndices.Select(o => o.spin), next.SpinOrbitalIndices.Select(o => o.spin))
                        &&
                        Enumerable.SequenceEqual(curr.SpinOrbitalIndices.Select(o => o.orbital), next.SpinOrbitalIndices.Select(o => o.orbital)))
                    {
                        curr.AddCoeff(next.coeff);
                        output[nElements - 1] = curr;
                    }
                    else
                    {
                        output.Add(next);
                        nElements++;
                    }
                }
            }
            return output;
        }

        #endregion

        #region Comparers


        /// <summary>
        /// Compares two equal-length sequences of integers by perform an element-wise
        /// comparison, starting from the first element.
        /// </summary>
        /// <param name="xArr">First sequence of integers.</param>
        /// <param name="yArr">Second sequence of integers.</param>
        /// <returns>
        /// Returns <c>1</c> if <paramref name="xArr"/> is greater than <paramref name="yArr"/>.
        /// Returns <c>-1</c> if <paramref name="xArr"/> is less than <paramref name="yArr"/>.
        /// Returns <c>0</c> if <paramref name="xArr"/> is equal to <paramref name="yArr"/>.
        /// </returns>
        /// <example>
        /// CompareIntArray(new Int64[] {5}, new Int64[] {7}) == -1;
        /// CompareIntArray(new Int64[] {5,7}, new Int64[] {5,6}) == 1;
        /// CompareIntArray(new Int64[] {2,1,3}, new Int64[] {2,2,3}) == -1;
        /// </example>
        public static int CompareIntArray(IEnumerable<Int64> xArr, IEnumerable<Int64> yArr)
        {
            foreach (var item in xArr.Zip(yArr, (x, y) => (x, y)))
            {
                if (item.y > item.x)
                {
                    return -1;
                }
                else if (item.y < item.x)
                {
                    return 1;
                }
            }
            return 0;
        }

        /// <summary>
        /// IComparer for two integers.
        /// </summary>
        public class IntIComparer : IComparer<Int64>
        {
            public int Compare(Int64 x, Int64 y) => System.Math.Sign(x - y);
        }


        public class IntArrayIComparer : IComparer<IEnumerable<Int64>>
        {
            public int Compare(IEnumerable<Int64> x, IEnumerable<Int64> y)
            {
                return CompareIntArray(x, y);
            }
        }

        /// <summary>
        /// IComparer for <c>(QArray<Int64>, QArray<Double>)</c>. This compares only
        /// the integer sequence, and ignores the double sequence.
        /// </summary>
        public class HTermIndexIComparer : IComparer<HTerm>
        {
            public int Compare(HTerm x, HTerm y) =>
                CompareIntArray(x.Item1, y.Item1);
        }

        /// <summary>
        /// Equality comparer for <see cref="HTerm"/>. This compares both
        /// the integer sequence and the double sequence.
        /// </summary>
        public class HTermArrayComparer : IEqualityComparer<HTerm>
        {
            public bool Equals(HTerm x, HTerm y)
            {
                if (x.Item1.SequenceEqual(y.Item1) && x.Item2.SequenceEqual(y.Item2))
                {
                    return true;
                }
                else
                {
                    return false;
                }
            }
            public int GetHashCode(HTerm x)
            {
                int h = 19;
                foreach (var i in x.Item1.Select(o => o.GetHashCode()))
                {
                    h = h * 31 + i;
                }
                foreach (var i in x.Item2.Select(o => o.GetHashCode()))
                {
                    h = h * 17 + i;
                }
                return h;
            }
        }

        public class IntArrayIEqualityComparer : IEqualityComparer<IEnumerable<Int64>>
        {
            public bool Equals(IEnumerable<Int64> x, IEnumerable<Int64> y)
            {
                return x.SequenceEqual(y);
            }
            public int GetHashCode(IEnumerable<Int64> x)
            {
                int h = 19;
                foreach (var i in x)
                {
                    h = h * 31 + ((int)i);
                }
                return h;
            }
        }
        #endregion

        


    }
}
