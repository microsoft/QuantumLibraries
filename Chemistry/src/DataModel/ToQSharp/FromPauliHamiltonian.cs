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
    public static class ToQSharp
    {
        
        public static HTerm FromPauliTerm(PauliTerm term, PauliTermValue value)
        {
            return new HTerm((new QArray<Int64>(term.QubitIndices.Select(o=>(Int64)o)), new QArray<double>(value.Value)));
        }

        public static (double, int, JWOptimizedHTerms) ToQSharpFormat(this PauliHamiltonian pauliHamiltonian)
        {
            double energyOffset = pauliHamiltonian.Terms[TermType.PauliTerm.Identity].Values.First().Value.First();
            var nSpinOrbitals = pauliHamiltonian.SystemIndices.Count() + 1;
            var hZ = pauliHamiltonian.Terms[TermType.PauliTerm.Z].Select(o => FromPauliTerm(o.Key, o.Value)).ToList();
            var hZZ = pauliHamiltonian.Terms[TermType.PauliTerm.ZZ].Select(o => FromPauliTerm(o.Key, o.Value)).ToList();
            var hPQ = pauliHamiltonian.Terms[TermType.PauliTerm.PQ].Select(o => FromPauliTerm(o.Key, o.Value));
            var hPQQR = pauliHamiltonian.Terms[TermType.PauliTerm.PQQR].Select(o => FromPauliTerm(o.Key, o.Value));
            var hv0123 = pauliHamiltonian.Terms[TermType.PauliTerm.v01234].Select(o => FromPauliTerm(o.Key, o.Value)).ToList();
                                 
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
        ///public QArray<Int64> statePrep = new QArray<Int64>();    
        
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
    

}