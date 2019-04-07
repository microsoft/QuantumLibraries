// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Chemistry
{
    using System;

    /// <summary>
    /// All Hamiltonian terms must implement this interface.
    /// </summary>
    /// <typeparam name="TermClassification">Index to categories of terms.</typeparam>
    public interface ITermIndex <TermClassification>
    {
        TermClassification GetTermType();
    }

    /// <summary>
    /// All Hamiltonian terms must implement this interface.
    /// </summary>
    /// <typeparam name="TermClassification">Index to categories of terms.</typeparam>
    public interface ITermValue<TermValue>
    {
        TermValue AddValue(TermValue addThis);
        TermValue Default();

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
        public enum OrbitalIntegral
        {
            
            Identity, OneBody, TwoBody
        }

        public enum Fermion
        {
            Identity = 0, PP = 1, PQ = 2, PQQP = 3, PQQR = 4, PQRS = 5
        }

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
    /// Spin index enumeration type.
    /// </summary>
    public enum Spin : byte
    {
        u = 0, d = 1
    };

    /// <summary>
    /// Wavefunction types
    /// </summary>
    public enum StateType
    {
        Default = 0, SingleConfigurational = 1, SparseMultiConfigurational = 2, UnitaryCoupledCluster = 3
    }

    /// <summary>
    /// Enum for raising or lowering operator.
    /// </summary>
    public enum RaisingLowering
    {
        u = 0, d = 1, identity
    }

    public struct DoubleCoeff : ITermValue<DoubleCoeff>
    {
        public double Value;

        public DoubleCoeff(double value)
        {
            Value = value;
        }

        public static implicit operator DoubleCoeff(double value)
        {
            return new DoubleCoeff(value);
        }
        public static implicit operator double(DoubleCoeff value)
        {
            return value.Value;
        }

        public DoubleCoeff Default()
        {
            return new DoubleCoeff(0.0);
        }

        public DoubleCoeff AddValue(DoubleCoeff addThis)
        {
            return new DoubleCoeff(Value + addThis.Value);
        }

        /// <summary>
        /// Computes the L_p norm of term.
        /// </summary>
        /// <param name="power">Selects type of norm.</param>
        /// <returns>L_p norm of term.</returns>
        public double Norm(double power)
        {
            return Math.Abs(Value);
        }

        /// <summary>
        /// Override for string representation of Double
        /// </summary>
        /// <returns>String representation of Double</returns>
        public override string ToString()
        {
            return Value.ToString();
        }


        #region Equality Testing

        public override bool Equals(object obj)
        {
            return (obj is DoubleCoeff x) ? Equals(x) : false;
        }

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

        public override int GetHashCode()
        {
            return Value.GetHashCode();
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

        public static bool operator !=(DoubleCoeff x, DoubleCoeff y)
        {
            return !(x == y);
        }

        public static DoubleCoeff operator +(DoubleCoeff x, DoubleCoeff y)
        {
            return new DoubleCoeff(x.Value + y.Value);
        }

        public static DoubleCoeff operator -(DoubleCoeff x, DoubleCoeff y)
        {
            return new DoubleCoeff(x.Value - y.Value);
        }

        public static DoubleCoeff operator *(DoubleCoeff x, DoubleCoeff y)
        {
            return new DoubleCoeff(x.Value * y.Value);
        }


        #endregion

        /// <summary>
        ///     Represents the possible formats that can be used to represent integral
        ///     data sets.
        /// </summary>
        public enum IntegralDataFormat
        {
            LiQuiD, Broombridge
        }
    }
}