// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using Microsoft.Quantum.Simulation.Core;

using System;
using System.Linq;
using System.Collections.Generic;

using Microsoft.Quantum.Chemistry.Pauli;

namespace Microsoft.Quantum.Chemistry.ToQSharp
{
    using HTermArray = QArray<HTerm>;
    using Microsoft.Quantum.Chemistry.JordanWigner;
    
    public class FromPauliHamiltonian
    {
        public PauliHamiltonian pauliHamiltonian;

        public static HTerm ToQSharp(PauliTerm term, PauliTermValue value)
        {
            return new HTerm((new QArray<Int64>(term.QubitIndices.Select(o=>(Int64)o)), new QArray<double>(value.Value)));
        }

        public static (double, int, JWOptimizedHTerms) ToQSharp(PauliHamiltonian pauliHamiltonian)
        {
            double energyOffset = pauliHamiltonian.terms[TermType.PauliTerm.Identity].Values.First().Value.First();
            var nSpinOrbitals = pauliHamiltonian.systemIndices.Count() + 1;
            var hZ = pauliHamiltonian.terms[TermType.PauliTerm.Z].Select(o => ToQSharp(o.Key, o.Value));
            var hZZ = pauliHamiltonian.terms[TermType.PauliTerm.Z].Select(o => ToQSharp(o.Key, o.Value));
            var hPQ = pauliHamiltonian.terms[TermType.PauliTerm.Z].Select(o => ToQSharp(o.Key, o.Value));
            var hPQQR = pauliHamiltonian.terms[TermType.PauliTerm.Z].Select(o => ToQSharp(o.Key, o.Value));
            var hv0123 = pauliHamiltonian.terms[TermType.PauliTerm.Z].Select(o => ToQSharp(o.Key, o.Value));

            var hPQandPQQR = hPQ.Concat(hPQQR);

            return (
                energyOffset,
                nSpinOrbitals,
                new JWOptimizedHTerms((
                new HTermArray(hZ),
                new HTermArray(hZZ),
                new HTermArray(hPQandPQQR),
                new HTermArray(hv0123)))
                );
        }
        /*
        /// <summary>
        /// Returns data for consumption by Q#.
        /// </summary>
        public static JordanWignerEncodingData ToQSharp(
            PauliHamiltonian pauliHamiltonian,
            string selectInputState = "Greedy")
        {
            var (energyOffset, nSpinOrbitals, hamiltonianData) = ToQSharp(pauliHamiltonian);

            return new JordanWignerEncodingData(
                (nSpinOrbitals
                , hamiltonianData
                , InputStateToQSharp(InputStates[selectInputState])
                , energyOffset));
        }
        */
        // Number of orbitals
        /// <summary>
        /// <see cref="NOrbitals"/> is the number of distinct orbitals.
        /// <see cref="NSpinOrbitals"/> is the number of distinct spin orbitals.
        /// </summary>
        public Int64 NOrbitals, NSpinOrbitals;
        /// <summary>
        /// This indexes the spin-orbitals to be occupied by electrons.
        /// </summary>
        public QArray<Int64> statePrep = new QArray<Int64>();    
        
        #region Term Sorting
        /// <summary>
        /// IComparer for sorting PQ and PQQR terms.
        /// </summary>
        public class HPQandPQQRIComparer : IComparer<HTerm>
        {
            public int Compare(HTerm x, HTerm y)
            {

                var xArr = x.Item1.Length == 2 ? new Int64[] { x.Item1[0], x.Item1[1], -1, -1 } : new Int64[] { x.Item1[0], x.Item1[3], x.Item1[2], x.Item1[1] };
                var yArr = y.Item1.Length == 2 ? new Int64[] { y.Item1[0], y.Item1[1], -1, -1 } : new Int64[] { y.Item1[0], y.Item1[3], y.Item1[2], y.Item1[1] };
                return Extensions.CompareIntArray(xArr, yArr);
            }
        }

        /// <summary>
        /// IComparer for <c>(QArray<Int64>, QArray<Double>)</c>. This compares only
        /// the integer sequence, and ignores the double sequence.
        /// </summary>
        public class HTermIndexIComparer : IComparer<HTerm>
        {
            public int Compare(HTerm x, HTerm y) =>
                Extensions.CompareIntArray(x.Item1, y.Item1);
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
        #endregion


    }
    
    public static class Extensions
    {


        #region Int64 Comparers


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
        public class Int64IComparer : IComparer<Int64>
        {
            public int Compare(Int64 x, Int64 y) => Math.Sign(x - y);
        }


        public class Int64ArrayIComparer : IComparer<IEnumerable<Int64>>
        {
            public int Compare(IEnumerable<Int64> x, IEnumerable<Int64> y)
            {
                return CompareIntArray(x, y);
            }
        }

        public class Int64ArrayIEqualityComparer : IEqualityComparer<IEnumerable<Int64>>
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