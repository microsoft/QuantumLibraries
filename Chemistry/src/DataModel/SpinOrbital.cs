// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.



using Microsoft.Quantum.Simulation.Core;

using System;
using System.Linq;
using System.Collections.Generic;

namespace Microsoft.Quantum.Chemistry
{
    /// <summary>
    /// Contains common data structures and methods relevant to representing a
    /// Fermion Hamiltonian.
    /// </summary>

    /// <summary>
    /// Spin index enumeration type.
    /// </summary>
    public enum Spin : Int64
    {
        u = 0, d = 1
    };

    /// <summary>
    /// Type representing a spin-orbital
    /// </summary>
    [Serializable]
    public struct SpinOrbital
    {
        /// <summary>
        /// Smallest orbital index
        /// </summary>
        internal const Int64 minOrbital = 0;

        /// <summary>
        /// Largest orbital index
        /// </summary>
        internal const Int64 maxOrbital = Int32.MaxValue - 1;

        /// <summary>
        /// Smallest spin index
        /// </summary>
        internal const Int64 minSpin = 0;

        /// <summary>
        /// Largest spin index
        /// </summary>
        internal const Int64 maxSpin = Int32.MaxValue - 1;

        /// <summary>
        /// <c>orbital</c> is the orbital index. This integer must satisfy '0 &lt;= orbital &lt Int32.ManValue'.
        /// </summary>
        public Int64 orbital;
        /// <summary>
        /// <c>spin</c> is the spin index. For electrons, spin=0 is spin up, and spin=1 is spin down. 
        /// This integer must satisfy '0 &lt;= orbital &lt Int16.MaxValue'.
        /// </summary>
        public Int64 spin;

        /// <summary>
        /// This maps the separate orbital and spin indices of a <c>SpinOrbital</c> 
        /// to a single integer index. 
        /// </summary>
        /// <param name="nOrbitals">The total number of orbitals.</param>
        public Int64 ToInt(Int64 nOrbitals = maxOrbital)
        {
           if (Settings.SetIndexConvention == Settings.IndexConventionType.UpDown)
            {
                return orbital * 2 + spin;
            }
            else
            {
                return orbital + nOrbitals * spin;
            }
        }


        /// <summary>
        /// This maps the single integer index produced by <c>ToInt(Int64 nOrbitals)</c>
        /// back to an orbital and spin index.
        /// </summary>
        /// <param name="nOrbitals">The total number of orbitals.</param>
        /// <param name="spinOrbitalIdx">A single integer representing a spin-orbital.</param>
        public SpinOrbital(Int64 nOrbitals, Int64 spinOrbitalIdx)
        {
            orbital = spinOrbitalIdx % nOrbitals;
            spin = spinOrbitalIdx / nOrbitals;
            IsValid();
        }

        /// <summary>
        /// Spin-orbital constructor.
        /// </summary>
        /// <param name="orbitalIdx">Orbital index.</param>
        /// <param name="spinIdx">Spin index.</param>
        public SpinOrbital(Int64 orbitalIdx, Spin spinIdx = Spin.u)
        {
            orbital = orbitalIdx;
            spin = (Int64)spinIdx;
            IsValid();
        }

        /// <summary>
        /// Spin-orbital constructor.
        /// </summary>
        /// <param name="idx">Tuple of (orbital index, spin index).</param>
        public SpinOrbital((Int64, Spin) idx)
        {
            orbital = idx.Item1;
            spin = (Int64)idx.Item2;
            IsValid();
        }

        /// <summary>
        /// Spin-orbital constructor.
        /// </summary>
        /// <param name="idx">Tuple of (orbital index, spin index).</param>
        public SpinOrbital((Int64, Int64) idx)
        {
            orbital = idx.Item1;
            spin = idx.Item2;
            IsValid();
        }

        /// <summary>
        /// Throws an exception if spin-orbital is invalid.
        /// </summary>
        /// <returns>Returns true if spin-orbital is valid.</returns>
        public bool IsValid()
        {
            if (orbital < minOrbital || orbital > maxOrbital || spin < minSpin || spin > maxSpin)
            {
                throw new System.ArgumentException(
                    $"Invalid SpinOrbital specified. " +
                    $"Orbital index is {orbital}. " +
                    $"It must satisfy {minOrbital} <= {orbital} <= {maxOrbital}." +
                    $"Spin index is {spin}. " +
                    $"It must satisfy {minSpin} <= {spin} <= {maxSpin}."
                    );
            }
            else
            {
                return true;
            }
        }

        /// <summary>
        /// This maps an array of integer indices to an array of SpinOrbitals.
        /// </summary>
        /// <param name="nOrbitals">The total number of orbitals.</param>
        /// <param name="spinOrbitalIdx">A sequence of integers representing a sequence of spin-orbitals.</param>
        public static SpinOrbital[] ToSpinOrbitals(Int64 nOrbitals, IEnumerable<Int64> fermionIdxArray)
        {
            return fermionIdxArray.Select(x => new SpinOrbital(nOrbitals, x)).ToArray();
        }

        /// <summary>
        /// This creates all possible spin-orbitals Hamiltonian terms given 
        /// an array of an even number of orbital indices. 
        /// </summary>
        /// <param name="nSpins">The number of possible spin states. This is 2 for electrons (spin 1/2).</param>
        /// <param name="orbitals">A sequence of integers representing a sequence of orbitals in Mullikan convention.</param>
        public static IEnumerable<SpinOrbital[]> Enumerate(Int64 nSpins, Int64[] orbitals)
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
                    spinOrbitalArray[idxOrbital] = new SpinOrbital
                    {
                        orbital = orbitals[idxOrbital],
                        spin = (idx / nSpins.Pow(idxOrbital)) % nSpins
                    };
                }
                yield return spinOrbitalArray;
            }
        }

        /// <summary>
        /// Override for string representation of spin-orbital data.
        /// </summary>
        /// <returns>String representation of spin-orbital data.</returns>
        public override string ToString()
        {
            return $"({ orbital},{ spin })";
        }

        /// <summary>
        /// Boolean equality operator definition.
        /// </summary>
        public static bool operator ==(SpinOrbital x, SpinOrbital y)
        {
            return x.orbital == y.orbital && x.spin == y.spin;
        }
        /// <summary>
        /// Boolean inequality operator definition.
        /// </summary>
        public static bool operator !=(SpinOrbital x, SpinOrbital y)
        {
            return !(x == y);
        }
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

        public override bool Equals(object x)
        {
            return Equals((SpinOrbital) x);
        }

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

        public override int GetHashCode()
        {
            return ToInt().GetHashCode();
        }
    }
}

