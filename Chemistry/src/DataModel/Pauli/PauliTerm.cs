// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using System;
using System.Linq;
using System.Collections.Generic;


namespace Microsoft.Quantum.Chemistry.Pauli
{
    public enum QubitEncoding
    {
        JordanWigner
    }
    
    /// <summary>
    /// Data structure for sparse representations of Pauli terms.
    /// </summary>
    public struct PauliTerm : 
        ITermIndex<TermType.PauliTerm>, 
        IEquatable<PauliTerm>
    {


        /// <summary>
        /// Qubit indices that represent this Pauli string.
        /// </summary>
        public List<int> QubitIndices { get; set; }

        /// <summary>
        /// LadderType of Pauli string encoded by list of qubit indices.
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
        /// Returns the sign of this term.
        /// </summary>
        /// <returns>Sign of this term.</returns>
        public int GetSign()
        {
            return 1;
        }

        /// <summary>
        /// Sets the sign of this fermion term to one.
        /// </summary>
        public void ResetSign()
        {
            // All Pauli terms currently have no sign.
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

    public struct PauliTTermValue : ITermValue<PauliTTermValue>
    {
        public double[] Value;

        public PauliTTermValue(double value)
        {
            Value = new[] { value };
        }

        public PauliTTermValue(IEnumerable<double> value)
        {
            Value = value.ToArray();
        }

        public PauliTTermValue Default()
        {
            return new PauliTTermValue(0.0);
        }

        public PauliTTermValue AddValue(PauliTTermValue addThis, int sign = 1)
        {
            return new PauliTTermValue( Value.Zip(addThis.Value, (a, b) => (a + b) * (double) sign));
        }

        public override string ToString()
        {
            return string.Join(", ", Value);
        }

        /// <summary>
        /// Computes the L_p norm of term.
        /// </summary>
        /// <param name="power">Selects type of norm.</param>
        /// <returns>L_p norm of term.</returns>
        public double Norm(double power)
        {
            return Value.Norm(power);
        }

        #region Equality Testing

        public override bool Equals(object obj)
        {
            return (obj is PauliTTermValue x) ? Equals(x) : false;
        }

        public bool Equals(PauliTTermValue x)
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
            return Value == x.Value;
        }

        public override int GetHashCode()
        {
            return Value.GetHashCode();
        }

        public static bool operator ==(PauliTTermValue x, PauliTTermValue y)
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

        public static bool operator !=(PauliTTermValue x, PauliTTermValue y)
        {
            return !(x == y);
        }

        public static PauliTTermValue operator +(PauliTTermValue x, PauliTTermValue y)
        {
            return new PauliTTermValue(x.Value.Zip(y.Value, (a, b) => (a + b)));
        }

        public static PauliTTermValue operator -(PauliTTermValue x, PauliTTermValue y)
        {
            return new PauliTTermValue(x.Value.Zip(y.Value, (a, b) => (a - b)));
        }

        public static PauliTTermValue operator *(PauliTTermValue x, PauliTTermValue y)
        {
            return new PauliTTermValue(x.Value.Zip(y.Value, (a, b) => (a * b)));
        }
        #endregion
    }
  
}

