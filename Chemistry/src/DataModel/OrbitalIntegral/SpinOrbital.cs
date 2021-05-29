// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using System;
using System.Collections.Generic;

namespace Microsoft.Quantum.Chemistry.OrbitalIntegrals
{
    /// <summary>
    /// Indexing scheme representing a spin-orbital.
    /// </summary>
    public class SpinOrbital : 
        IEquatable<SpinOrbital>, 
        IComparable<SpinOrbital>
    {

        IndexConvention IdxConvention;

        /// <summary>
        /// Smallest orbital index
        /// </summary>
        internal const int minOrbital = 0;

        /// <summary>
        /// Largest orbital index
        /// </summary>
        internal const int maxOrbital = int.MaxValue - 1;

        /// <summary>
        /// Smallest spin index
        /// </summary>
        internal const int minSpin = 0;

        /// <summary>
        /// Largest spin index
        /// </summary>
        internal const int maxSpin = int.MaxValue - 1;

        /// <summary>
        /// <c>orbital</c> is the orbital index. This integer must satisfy '0 &lt;= orbital &lt; Int32.ManValue'.
        /// </summary>
        public int Orbital;
        /// <summary>
        /// <c>spin</c> is the spin index. For electrons, spin=0 is spin up, and spin=1 is spin down. 
        /// This integer must satisfy '0 &lt;= orbital &lt; Int16.MaxValue'.
        /// </summary>
        public int Spin;

        /// <summary>
        /// This maps the separate orbital and spin indices of a <c>SpinOrbital</c> 
        /// to a single integer index. 
        /// </summary>
        /// <param name="nOrbitals">The total number of orbitals.</param>
        public int ToInt(IndexConvention indexConvention, int nOrbitals = maxOrbital) =>
            indexConvention == IndexConvention.UpDown ? Orbital * 2 + Spin : Orbital + nOrbitals * Spin;

        /// <summary>
        /// This maps the separate orbital and spin indices of a <c>SpinOrbital</c> 
        /// to a single integer index according to `2*orbitalIndex + spinIndex`.
        /// </summary>
        public int ToInt() => ToInt(IndexConvention.UpDown);

        // This exists as a convenience function for creating spin-orbitals in samples.
        /// <summary>
        /// Implicit operator for creating a spin-orbital.
        /// </summary>
        /// <param name="setIndex">Tuple where the first parameter
        /// is the orbital index, and the second parameter
        /// is the spin index.</param>
        public static implicit operator SpinOrbital((int, Spin) setIndex)
        {
            return new SpinOrbital(setIndex);
        }

        /// <summary>
        /// Spin-orbital constructor.
        /// </summary>
        /// <param name="orbitalIdx">Orbital index.</param>
        /// <param name="spinIdx">Spin index.</param>
        public SpinOrbital(int orbitalIdx, Spin spinIdx)
        {
            Orbital = orbitalIdx;
            Spin = (int)spinIdx;
            ThrowIfInvalid();
        }

        /// <summary>
        /// Spin-orbital constructor.
        /// </summary>
        /// <param name="idx">Tuple of (orbital index, spin index).</param>
        public SpinOrbital((int, Spin) idx) : this(idx.Item1, idx.Item2) { }
            
        /// <summary>
        /// Spin-orbital constructor.
        /// </summary>
        /// <param name="idx">Tuple of (orbital index, spin index).</param>
        public SpinOrbital((int, int) idx) : this(idx.Item1, (Spin) idx.Item2) { }

        /// <summary>
        /// Empty spin-orbital constructor.
        /// </summary>
        public SpinOrbital()
        {

        }

        /// <summary>
        /// Throws an exception if spin-orbital is invalid.
        /// </summary>
        /// <returns>Returns true if spin-orbital is valid.</returns>
        public void ThrowIfInvalid()
        {
            if (Orbital < minOrbital || Orbital > maxOrbital || Spin < minSpin || Spin > maxSpin)
            {
                throw new System.ArgumentException(
                    $"Invalid SpinOrbital specified. " +
                    $"Orbital index is {Orbital}. " +
                    $"It must satisfy {minOrbital} <= {Orbital} <= {maxOrbital}." +
                    $"Spin index is {Spin}. " +
                    $"It must satisfy {minSpin} <= {Spin} <= {maxSpin}."
                    );
            }
        }

        /// <summary>
        /// This creates all possible spin-orbitals Hamiltonian terms given 
        /// an array of an even number of orbital indices. 
        /// </summary>
        /// <param name="nSpins">The number of possible spin states. This is 2 for electrons (spin 1/2).</param>
        /// <param name="orbitals">A sequence of integers representing a sequence of orbitals in Mullikan convention.</param>
        public static IEnumerable<SpinOrbital[]> Enumerate(int[] orbitals, int nSpins = 2)
        {
            // Assumes spinOrbitals has an even number of elements
            // Only index over like spins S1 S2 S3 ... S3 S2 S1
            var orbitalArrayLength = orbitals.Length;
            var nSpinOrbitalArrays = nSpins.Pow(orbitalArrayLength);

            for (int idx = 0; idx < nSpinOrbitalArrays; idx++)
            {
                SpinOrbital[] spinOrbitalArray = new SpinOrbital[orbitalArrayLength];
                for (int idxOrbital = 0; idxOrbital < orbitalArrayLength; idxOrbital++)
                {
                    spinOrbitalArray[idxOrbital] = new SpinOrbital(orbitals[idxOrbital], (Spin)((idx / nSpins.Pow(idxOrbital)) % nSpins));
                }
                yield return spinOrbitalArray;
            }
        }
        
        /// <summary>
        /// Override for string representation of spin-orbital data.
        /// </summary>
        /// <returns>String representation of spin-orbital data.</returns>
        public override string ToString() => $"({ Orbital},{ Spin })";

        /// <summary>
        /// Boolean equality operator definition.
        /// </summary>
        public static bool operator == (SpinOrbital x, SpinOrbital y) => x.Orbital == y.Orbital && x.Spin == y.Spin;

        /// <summary>
        /// Boolean inequality operator definition.
        /// </summary>
        public static bool operator != (SpinOrbital x, SpinOrbital y) => !(x == y);

        /// <summary>
        /// Default comparer for spin orbitals first compares the orbital
        /// index, then second compares the spin index.
        /// </summary>
        /// <param name="x">Spin orbital to compare with.</param>
        /// <returns>Result of the comparison.</returns>
        public int CompareTo(SpinOrbital x) => ToInt().CompareTo(x.ToInt());

        /*
        /// <summary>
        /// Boolean greater-than operator definition.
        /// </summary>
        public static bool operator >(SpinOrbital x, SpinOrbital y)
        {
            return x.ToInt() > y.ToInt();
        }
        /// <summary>
        /// Boolean less-than operator definition.
        /// </summary>
        public static bool operator <(SpinOrbital x, SpinOrbital y)
        {
            return x.ToInt() < y.ToInt();
        }
        */
        public override bool Equals(object x) => Equals((SpinOrbital) x);

        public bool Equals(SpinOrbital x)
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

        public override int GetHashCode() => ToInt(IndexConvention.UpDown).GetHashCode();
    }
}

