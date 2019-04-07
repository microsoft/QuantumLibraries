// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using Microsoft.Quantum.Simulation.Core;

using System;
using System.Linq;
using System.Collections.Generic;

namespace Microsoft.Quantum.Chemistry
{
   
    
    /// <summary>
    /// Data structure for sparse representations of Pauli terms.
    /// </summary>
    public struct PauliTerm : 
        HamiltonianTerm<TermType.PauliTerm>, 
        IEquatable<PauliTerm>
    {
        public enum Encoding
        {
            JordanWigner
        }

        /// <summary>
        /// Qubit indices that represent this Pauli string.
        /// </summary>
        public List<int> QubitIndices { get; set; }

        /// <summary>
        /// Type of Pauli string encoded by list of qubit indices.
        /// </summary>
        public TermType.PauliTerm TermType;

        #region Constructors
        /// <summary>
        /// Construct a copy of the input instance.
        /// </summary>
        /// <param name="pauliString">Input.</param>
        public PauliTerm(PauliTerm pauliString)
        {
            // All constructions are pass by value.
            QubitIndices = pauliString.QubitIndices.ToList();
            TermType = pauliString.TermType;
        }

        /// <summary>
        /// Construct instance from sequence of qubit indices.
        /// </summary>
        /// <param name="pauliString">Input.</param>
        public PauliTerm(IEnumerable<int> qubitIndices, TermType.PauliTerm type)
        {
            QubitIndices = qubitIndices.ToList();
            TermType = type;
        }
        #endregion
        
        /// <summary>
        /// Return the category of this fermion term.
        /// </summary>
        /// <returns>Category of fermion term.</returns>
        public TermType.PauliTerm GetTermType()
        {
            return TermType;
        }

        /// <summary>
        /// Returns a human-readable description this object.
        /// </summary>
        public override string ToString()
        {
            return $"{TermType}: [ {string.Join(" ", QubitIndices)} ]";
        }

        #region Equality Testing

        public override bool Equals(object obj)
        {
            return (obj is PauliTerm x) ? Equals(x) : false;
        }

        public bool Equals(PauliTerm x)
        {
            // If parameter is null, return false.
            if (ReferenceEquals(x, null))
            {
                return false;
            }

            // Optimization for a common success case.
            if (ReferenceEquals(this, x))
            {
                return true;
            }

            // If run-time types are not exactly the same, return false.
            if (GetType() != x.GetType())
            {
                return false;
            }
            // Return true if the fields match.
            return QubitIndices.SequenceEqual(x.QubitIndices) && TermType == x.TermType;
        }

        public override int GetHashCode()
        {
            int h = 19;
            foreach (var i in QubitIndices)
            {
                h = h * 53 + i.GetHashCode() + 97 * (int) TermType;
            }
            return h;
        }

        public static bool operator == (PauliTerm x, PauliTerm y)
        {
            // Check for null on left side.
            if (Object.ReferenceEquals(x, null))
            {
                if (Object.ReferenceEquals(y, null))
                {
                    // null == null = true.
                    return true;
                }

                // Only the left side is null.
                return false;
            }
            // Equals handles case of null on right side.
            return x.Equals(y);
        }

        public static bool operator !=(PauliTerm x, PauliTerm y)
        {
            return !(x == y);
        }
        #endregion

    }
}

