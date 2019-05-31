// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using System;
using System.Linq;
using System.Collections.Generic;


namespace Microsoft.Quantum.Chemistry.Paulis
{
    public enum QubitEncoding
    {
        JordanWigner
    }
    
    /// <summary>
    /// Data structure for sparse representations of Pauli terms.
    /// </summary>
    public struct PauliTerm : ITermIndex<TermType.PauliTerm, PauliTerm> 
    {


        /// <summary>
        /// Qubit indices that represent this Pauli string.
        /// </summary>
        public int[] QubitIndices { get; set; }

        /// <summary>
        /// LadderType of Pauli string encoded by list of qubit indices.
        /// </summary>
        public TermType.PauliTerm TermType { get; set; }

        #region Constructors
        /// <summary>
        /// Construct a copy of the input instance.
        /// </summary>
        /// <param name="pauliString">Input.</param>
        public PauliTerm(PauliTerm pauliString)
        {
            // All constructions are pass by value.
            QubitIndices = pauliString.QubitIndices.Clone<int>();
            TermType = pauliString.TermType;
        }

        /// <summary>
        /// Construct instance from sequence of qubit indices.
        /// </summary>
        /// <param name="pauliString">Input.</param>
        public PauliTerm(IEnumerable<int> qubitIndices, TermType.PauliTerm type)
        {
            QubitIndices = qubitIndices.ToArray().Clone<int>();
            TermType = type;
        }
        #endregion
        

        /// <summary>
        /// Returns the sign of this term.
        /// </summary>
        /// <returns>Sign of this term.</returns>
        public int Sign => 1;

        /// <summary>
        /// Sets the sign of this fermion term to one.
        /// </summary>
        public void ResetSign()
        {
            // All Pauli terms currently have no sign.
        }

        /// <summary>
        /// Returns a human-readable description of this object.
        /// </summary>
        public override string ToString() => $"{TermType}: [ {string.Join(" ", QubitIndices)} ]";

        /// <summary>
        /// Creates a copy of this instance.
        /// </summary>
        public PauliTerm Clone() => new PauliTerm(QubitIndices.ToArray().Clone<int>(), TermType);

        #region Equality Testing

        public override bool Equals(object obj) => (obj is PauliTerm x) ? Equals(x) : false;

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

        public static bool operator !=(PauliTerm x, PauliTerm y) => !(x == y);
        #endregion

    }

    public struct PauliTermValue : ITermValue<PauliTermValue>
    {
        public double[] Value;

        public PauliTermValue(double value)
        {
            Value = new[] { value };
        }

        public PauliTermValue(IEnumerable<double> value)
        {
            Value = value.ToArray().Clone<double>();
        }

        public PauliTermValue Default() => new PauliTermValue(0.0);

        public PauliTermValue AddValue(PauliTermValue addThis, int sign = 1) =>
            new PauliTermValue( Value.Zip(addThis.Value, (a, b) => a + (b * (double) sign)));
        
        public PauliTermValue SetValue(PauliTermValue setThis, int sign = 1) =>
            new PauliTermValue(setThis.Value.Select( a => a * (double)sign));

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

        /// <summary>
        /// Creates a copy of this instance.
        /// </summary>
        public PauliTermValue Clone() => new PauliTermValue(Value.Clone<double>());

        #region Equality Testing

        public override bool Equals(object obj)
        {
            return (obj is PauliTermValue x) ? Equals(x) : false;
        }

        public bool Equals(PauliTermValue x)
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

        public static bool operator ==(PauliTermValue x, PauliTermValue y)
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

        public static bool operator !=(PauliTermValue x, PauliTermValue y)
        {
            return !(x == y);
        }

        public static PauliTermValue operator +(PauliTermValue x, PauliTermValue y)
        {
            return new PauliTermValue(x.Value.Zip(y.Value, (a, b) => (a + b)));
        }

        public static PauliTermValue operator -(PauliTermValue x, PauliTermValue y)
        {
            return new PauliTermValue(x.Value.Zip(y.Value, (a, b) => (a - b)));
        }

        public static PauliTermValue operator *(PauliTermValue x, PauliTermValue y)
        {
            return new PauliTermValue(x.Value.Zip(y.Value, (a, b) => (a * b)));
        }
        #endregion
    }
  
}

