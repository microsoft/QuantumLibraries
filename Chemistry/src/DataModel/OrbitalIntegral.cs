// Copyright (c) Microsoft Corporation. All rights reserved. Licensed under the
// Microsoft Software License Terms for Microsoft Quantum Simulation Library (Preview).
// See LICENSE.md in the project root for license information.


using Microsoft.Quantum.Simulation.Core;

using System;
using System.Linq;
using System.Collections.Generic;

namespace Microsoft.Quantum.Chemistry
{


    /// <summary>
    /// Type representing orbital overlap integrals.
    /// </summary>
    public struct OrbitalIntegral
    {
        public enum Convention
        {
            Dirac, Mulliken
        }

        /// <summary>
        /// <c>Int64[] OrbitalIndices</c> represents the indices of orbitals in the overlap integral.
        /// </summary>
        public Int64[] OrbitalIndices;

        /// <summary>
        /// <c>Double coefficient</c> represents the coefficient of the orbital overlap integral.
        /// </summary>
        public Double Coefficient;

        /// <summary>
        /// OrbitalIntegral constructor.
        /// </summary>
        /// <param name="orbitalIndices">Array of orbital indices in Dirac notation.</param>
        /// <param name="coefficient">coefficient of orbital integral.</param>
        public OrbitalIntegral(IEnumerable<Int64> orbitalIndices, Double coefficient = 0.0) 
        {
            OrbitalIndices = orbitalIndices.ToArray();
            Coefficient = coefficient;
        }

        /// <summary>
        /// OrbitalIntegral constructor.
        /// </summary>
        /// <param name="orbitalIndices">Array of orbital indices in Dirac notation.</param>
        /// <param name="coefficient">coefficient of orbital integral.</param>
        public OrbitalIntegral(IEnumerable<int> orbitalIndices, Double coefficient = 0.0) :
            this(orbitalIndices.Select(o => (Int64)o).ToArray(), coefficient)
        { }

        /// <summary>
        /// OrbitalIntegral constructor.
        /// </summary>
        /// <param name="orbitalIndices">Array of orbital indices.</param>
        /// <param name="coefficient">coefficient of orbital integral.</param>
        /// <param name="convention">Convention for ordering of orbital indices.</param>
        public OrbitalIntegral(IEnumerable<Int64> orbitalIndices, Double coefficient, Convention convention)
        {
            if (convention == Convention.Mulliken)
            {
                if (orbitalIndices.Count() == 2)
                {
                    OrbitalIndices = orbitalIndices.Select(o => (Int64)o).ToArray();
                }
                else if (orbitalIndices.Count() == 4)
                {
                    var p = orbitalIndices.ElementAt(0);
                    var q = orbitalIndices.ElementAt(2);
                    var r = orbitalIndices.ElementAt(3);
                    var s = orbitalIndices.ElementAt(1);
                    OrbitalIndices = new Int64[] { p, q, r, s };
                }
                else
                {
                    throw new System.ArgumentException("Mulliken convention for 2 or 4 indices is not defined.");
                }
                Coefficient = coefficient;
            }
            else
            {
                OrbitalIndices = orbitalIndices.ToArray();
                Coefficient = coefficient;
            }
        }

        /// <summary>
        /// Returns length of <see cref="OrbitalIndices"/>. 
        /// </summary>
        /// <returns>Length of <see cref="OrbitalIndices"/>.</returns>
        public Int64 Length()
        {
            return OrbitalIndices.Length;
        }


        /// <summary>
        /// Enumerates over all orbital integrals with the same coefficient
        /// as this.
        /// These symmetries arise from: 
        /// - Indistinguishability of electrons.
        /// - Orbitals are assumed to be real.
        /// </summary>
        /// <returns>Array of orbital integrals with the same coefficient.</returns>
        public OrbitalIntegral[] EnumerateOrbitalSymmetries()
        {
            var coefficient = Coefficient;
            if (OrbitalIndices.Length == 2)
            {
                var i = OrbitalIndices[0];
                var j = OrbitalIndices[1];
                var symmetries = new Int64[][]
                {
                    new Int64[] {i, j},
                    new Int64[] {j, i}
                };
                return symmetries.Distinct(new Extensions.IntArrayIEqualityComparer()).Select(o => new OrbitalIntegral(o, coefficient)).ToArray();
            }
            else if (OrbitalIndices.Length == 4)
            {
                var i = OrbitalIndices[0];
                var j = OrbitalIndices[1];
                var k = OrbitalIndices[2];
                var l = OrbitalIndices[3];
                var symmetries = new Int64[][] {
                    new Int64[] { i, j, k, l }, // 0123
                    new Int64[] { j, i, l, k }, // 1032
                    new Int64[] { k, l, i, j }, // 2301
                    new Int64[] { l, k, j, i }, // 3210
                    new Int64[] { i, k, j, l }, // 0213
                    new Int64[] { k, i, l, j }, // 2031
                    new Int64[] { j, l, i, k }, // 1302
                    new Int64[] { l, j, k, i }  // 3120
                };
                return symmetries.Distinct(new Extensions.IntArrayIEqualityComparer()).Select(o => new OrbitalIntegral(o, coefficient)).ToArray();
            }
            else
            {
                throw new System.NotImplementedException();
            }
        }

        /// <summary>
        /// Returns orbital indices sorted in a canonical form that generates
        /// the same set of orbital integrals through <see cref="EnumerateSpinOrbitals"/>.
        /// </summary>
        /// <returns>An <see cref="OrbitalIntegral"/> in canonical form.</returns>
        public OrbitalIntegral ToCanonicalForm()
        {
            var symmetries = EnumerateOrbitalSymmetries().Select(o => o.OrbitalIndices).ToList();
            symmetries.Sort(new Extensions.IntArrayIComparer());
            return new OrbitalIntegral(symmetries.First(), Coefficient);
        }

        /// <summary>
        /// Checks of this <see cref="OrbitalIntegral"/> has indices sorted in canonical order.
        /// </summary>
        /// <returns>Returns <see cref="bool"/> if <see cref="OrbitalIntegral"/> is canonically sorted
        /// and <see cref="false"/> otherwise.
        /// </returns>
        public bool IsInCanonicalOrder()
        {
            if (Length() == 2 || Length() == 4)
            {
                var canonicalOrder = ToCanonicalForm();
                return canonicalOrder.OrbitalIndices.SequenceEqual(OrbitalIndices);
            }
            else
            {
                throw new System.NotImplementedException();
            }
        }

        /// <summary>
        /// Override for string representation of orbital integral data.
        /// </summary>
        /// <returns>String representation of orbital integral data.</returns>
        public override string ToString()
        {
            return $"([{String.Join(",", OrbitalIndices)}], {Coefficient})";
        }

        /// <summary>
        /// Enumerates over all spin-orbital integrals represented 
        /// by this term.
        /// These symmetries arise from: 
        /// - Indistinguishability of electrons.
        /// - Orbitals are assumed to be real.
        /// - electron spins must be paired.
        /// </summary>
        /// <returns>
        /// Array of spin-orbital orbital integrals with the same coefficient.
        /// </returns>
        public SpinOrbital[][] EnumerateSpinOrbitals()
        {
            // Assumes spinOrbitals has an even number of elements
            // Only index over like spins S1 S2 S3 ... S3 S2 S1
            const Int64 nSpins = 2L;
            var nSpinOrbitalArrays = nSpins.Pow(OrbitalIndices.Length / 2);
            SpinOrbital[][] spinOrbitalArrayOfArray = new SpinOrbital[OrbitalIndices.Length][];
            for (int idx = 0; idx < nSpinOrbitalArrays; idx++)
            {
                SpinOrbital[] spinOrbitalArray = new SpinOrbital[OrbitalIndices.Length];
                for (int idxOrbital = 0; idxOrbital < OrbitalIndices.Length / 2; idxOrbital++)
                {
                    var fst = idxOrbital;
                    var lst = OrbitalIndices.Length - idxOrbital - 1;
                    var spin = (idx / nSpins.Pow(idxOrbital)) % nSpins;
                    spinOrbitalArray[fst] = new SpinOrbital
                    {
                        orbital = OrbitalIndices[fst],
                        spin = spin
                    };
                    spinOrbitalArray[lst] = new SpinOrbital
                    {
                        orbital = OrbitalIndices[lst],
                        spin = spin
                    };
                }
                spinOrbitalArrayOfArray[idx] = spinOrbitalArray;
            }
            return spinOrbitalArrayOfArray;
        }
    }
}

