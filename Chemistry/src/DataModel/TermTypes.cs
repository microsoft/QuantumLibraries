// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Chemistry
{
    using System;

    using Newtonsoft.Json;
    using Newtonsoft.Json.Converters;

    /// <summary>
    /// All Hamiltonian terms must implement this interface.
    /// </summary>
    /// <typeparam name="TTermClassification">Index to categories of terms.</typeparam>
    public interface ITermIndex <TTermClassification, TTermIndex> : IEquatable<TTermIndex>
    {
        /// <summary>
        /// Compute the classification of the Hamiltonian term.
        /// </summary>
        /// <returns>Classification of Hamiltonian term.</returns>
        TTermClassification TermType { get; }

        /// <summary>
        /// Clones the values of this object.
        /// </summary>
        /// <returns>A deep copy of this object.</returns>
        TTermIndex Clone();
        /// <summary>
        /// Obtain the sign of the Hamiltonian term, if any.
        /// </summary>
        /// <returns>Sign of the term.</returns>
        int Sign { get; }

        /// <summary>
        /// Sets the sign of the Hamiltonian term to one.
        /// </summary>
        void ResetSign();
    }

    /// <summary>
    /// All Hamiltonian terms must implement this interface.
    /// </summary>
    /// <typeparam name="TTermClassification">Index to categories of terms.</typeparam>
    public interface ITermValue<TTermValue> : IEquatable<TTermValue>
    {
        /// <summary>
        /// Sets the value of the Hamiltonian term.
        /// </summary>
        /// <param name="setThis">Desired value of term.</param>
        /// <param name="sign">Multiply the applied value by this parameter.</param>
        /// <returns></returns>
        TTermValue SetValue(TTermValue setThis, int sign);
        TTermValue AddValue(TTermValue addThis, int sign);
        TTermValue Default();

        TTermValue Clone();


        /// <summary>
        /// Computes the L_p norm of term.
        /// </summary>
        /// <param name="power">Selects type of norm.</param>
        /// <returns>L_p norm of term.</returns>
        double Norm(double power);
    }



    /// <summary>
    /// Class containing a indices to a variety of term categories.
    /// </summary>
    // Might want to distribute these to the separate term classes.
    public static class TermType
    {
        [JsonConverter(typeof(StringEnumConverter))]
        public enum OrbitalIntegral
        {

            Identity, OneBody, TwoBody
        }

        [JsonConverter(typeof(StringEnumConverter))]
        public enum Fermion
        {
            Identity = 0, PP = 1, PQ = 2, PQQP = 3, PQQR = 4, PQRS = 5
        }

        [JsonConverter(typeof(StringEnumConverter))]
        public enum PauliTerm
        {
            Identity = 0,   // This has contribution from PP and PQQP terms.
            Z = 1,          // This has contribution from PP and PQQP terms.
            ZZ = 2,         // This has contribution from PQQP terms.
            PQ = 3,         // This has contributions from PQ and PQQR terms.
            PQQR = 4,       // This has contributions from PQQR terms.
            v01234 = 5      // This has contributions from PQRS terms.
        }
    }

    /// <summary>
    /// Spin index up/dpwn enumeration type.
    /// </summary>
    [JsonConverter(typeof(StringEnumConverter))]
    public enum Spin : byte
    {
        u = 0, d = 1
    };

    /// <summary>
    /// Available indexing convention from a spin-orbital index to an integer.
    /// </summary>
    public enum IndexConvention
    {
        NotApplicable, UpDown, HalfUp
    }

    /// <summary>
    /// Wavefunction types
    /// </summary>
    [JsonConverter(typeof(StringEnumConverter))]
    public enum StateType
    {
        NotRecognized = -1, Default = 0, SingleConfigurational = 1, SparseMultiConfigurational = 2, UnitaryCoupledCluster = 3
    }

    /// <summary>
    /// Enum for raising or lowering operator.
    /// </summary>
    [JsonConverter(typeof(StringEnumConverter))]
    public enum RaisingLowering
    {
        u = 0, d = 1, identity
    }

    /// <summary>
    ///     Represents the possible formats that can be used to represent integral
    ///     data sets.
    /// </summary>
    public enum IntegralDataFormat
    {
        /// <summary>
        /// Liquid format
        /// </summary>
        LiQuiD,
        /// <summary>
        /// Broombridge format
        /// </summary>
        Broombridge
    }


    /// <summary>
    /// Boxed version of double that implements the <see cref="ITermValue"></see> interface
    /// that requires all values to have a method to compute its norm.
    /// </summary>
    [JsonConverter(typeof(DoubleCoeff.JsonConverter))]
    public struct DoubleCoeff : ITermValue<DoubleCoeff>, IComparable
    {
        public double Value;

        public DoubleCoeff(double value)
        {
            Value = value;
        }

        public static implicit operator DoubleCoeff(double value) => new DoubleCoeff(value);
        public static implicit operator double(DoubleCoeff value) => value.Value;

        public DoubleCoeff Default() => 0.0;

        public DoubleCoeff AddValue(DoubleCoeff addThis, int sign) => Value + (addThis.Value * (double)sign);
        public DoubleCoeff SetValue(DoubleCoeff setThis, int sign) => setThis * (double)sign;

        /// <summary>
        /// Computes the L_p norm of term.
        /// </summary>
        /// <param name="power">Selects type of norm.</param>
        /// <returns>L_p norm of term.</returns>
        public double Norm(double power) => Math.Abs(Value);

        /// <summary>
        /// Override for string representation of Double
        /// </summary>
        /// <returns>String representation of Double</returns>
        public override string ToString() => Value.ToString();

        /// <summary>
        /// Creates a copy of this instance.
        /// </summary>
        public DoubleCoeff Clone()
        {
            return Value;
        }

        #region Equality Testing

        public override bool Equals(object obj) => (obj is DoubleCoeff x) ? Equals(x) : false;

        public bool Equals(DoubleCoeff x)
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

        public override int GetHashCode() => Value.GetHashCode();

        public int CompareTo(object obj)
        {
            if (obj is DoubleCoeff d) return this.Value.CompareTo(d.Value);

            return 0;
        }

        public static bool operator == (DoubleCoeff x, DoubleCoeff y)
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

        public static bool operator !=(DoubleCoeff x, DoubleCoeff y) => !(x == y);


        #endregion

        /// <summary>
        /// This JsonConverter encodes the DoubleCoeff as a double.
        /// </summary>
        public class JsonConverter : JsonConverter<DoubleCoeff>
        {
            /// <summary>
            /// Writers the LadderOperator as a (Type, Index) tuple.
            /// </summary>
            public override void WriteJson(JsonWriter writer, DoubleCoeff value, JsonSerializer serializer)
            {
                serializer.Serialize(writer, value.Value);
            }

            /// <summary>
            /// Reads the LadderOperator from a (Type, Index) tuple.
            /// </summary>
            public override DoubleCoeff ReadJson(JsonReader reader, Type objectType, DoubleCoeff existingValue, bool hasExistingValue, JsonSerializer serializer)
            {
                var value = serializer.Deserialize<Double>(reader);
                return new DoubleCoeff(value);
            }
        }

    }
}