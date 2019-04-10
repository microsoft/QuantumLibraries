// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using Microsoft.Quantum.Simulation.Core;

using System;
using System.Linq;
using System.Collections.Generic;

using Microsoft.Quantum.Chemistry.Pauli;
using Microsoft.Quantum.Chemistry.JordanWigner;

namespace Microsoft.Quantum.Chemistry.QSharpFormat
{
    /// <summary>
    /// Methods for converting electronic structure problem to data for consumption by Q#.
    /// </summary>
    public static partial class Convert
    {
        
        public static HTerm FromPauliTerm(PauliTerm term, PauliTTermValue value)
        {
            return new HTerm((new QArray<Int64>(term.QubitIndices.Select(o=>(Int64)o)), new QArray<double>(value.Value)));
        }

        internal static List<HTerm> CreateHTermList(PauliHamiltonian pauliHamiltonian, TermType.PauliTerm term) 
        {
            if (pauliHamiltonian.Terms.ContainsKey(term))
            {
                return pauliHamiltonian.Terms[term].Select(o => FromPauliTerm(o.Key, o.Value)).ToList();
            }
            else
            {
                return new List<HTerm>(); 
            }
        }

        public static (Double, Int64, JWOptimizedHTerms) ToQSharpFormat(this PauliHamiltonian pauliHamiltonian)
        {
            List<HTerm> DefaultHTerm = new List<HTerm>();

            var energyOffset = 0.0;
            if (pauliHamiltonian.Terms.ContainsKey(TermType.PauliTerm.Identity))
            {
                energyOffset = pauliHamiltonian.Terms[TermType.PauliTerm.Identity].Values.First().Value.First();
            }
            
            var nSpinOrbitals = pauliHamiltonian.SystemIndices.Max() + 1;
            var hZ = CreateHTermList(pauliHamiltonian,TermType.PauliTerm.Z);
            var hZZ = CreateHTermList(pauliHamiltonian,TermType.PauliTerm.ZZ);
            var hPQ = CreateHTermList(pauliHamiltonian,TermType.PauliTerm.PQ);
            var hPQQR = CreateHTermList(pauliHamiltonian,TermType.PauliTerm.PQQR);
            var hv0123 = CreateHTermList(pauliHamiltonian,TermType.PauliTerm.v01234);
                     

            var hPQandPQQR = hPQ.Concat(hPQQR).ToList();

            // Sort terms into desired order.
            hZ.Sort(new HTermIndexIComparer());
            hZZ.Sort(new HTermIndexIComparer());
            hPQandPQQR.Sort(new HPQandPQQRIComparer());
            hv0123.Sort(new HTermIndexIComparer());

            return (
                energyOffset,
                nSpinOrbitals,
                new JWOptimizedHTerms((
                new QArray<HTerm>(hZ),
                new QArray<HTerm>(hZZ),
                new QArray<HTerm>(hPQandPQQR),
                new QArray<HTerm>(hv0123)))
                );
        }

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
                return Extensions.CompareArray(xArr, yArr);
            }
        }

        /// <summary>
        /// IComparer for <c>(QArray<Int64>, QArray<Double>)</c>. This compares only
        /// the integer sequence, and ignores the double sequence.
        /// </summary>
        public class HTermIndexIComparer : IComparer<HTerm>
        {
            public int Compare(HTerm x, HTerm y) =>
                Extensions.CompareArray(x.Item1, y.Item1);
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
    

}