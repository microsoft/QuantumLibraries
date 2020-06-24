// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

#nullable enable

using System.Collections.Generic;
using Microsoft.Quantum.Chemistry.Broombridge;
using Microsoft.Quantum.Chemistry.Fermion;
using Microsoft.Quantum.Chemistry.OrbitalIntegrals;

namespace Microsoft.Quantum.Chemistry
{

    public static class QuantityExtensions
    {
        public static Quantity<TValue> WithUnits<TValue>(this TValue value, string units) =>
            new Quantity<TValue>
            {
                Units = units,
                Value = value
            };
    }

    public interface IHasUnits
    {
        public string Units { get; set; }
    }

    public interface IHasMetadata
    {
        public Dictionary<string, object> Metadata { get; set; }
    }

    public struct Quantity<TValue> : IHasUnits
    {
        public string Units { get ; set; }
        public TValue Value { get; set; }
    }

    public struct BoundedQuantity<TValue> : IHasUnits
    where TValue : struct
    {
        public string Units { get ; set; }
        public TValue? Value { get; set; }
        public TValue Lower { get; set; }
        public TValue Upper { get; set; }
    }

    public struct BasisSet
    {
        public string Type { get; set; }
        public string Name { get; set; }
    }

    
    public struct Geometry : IHasUnits
    {
        public string Units { get ; set; }
        public string CoordinateSystem { get; set; }

        public string Symmetry { get; set; }

        public List<Dictionary<string, object>> Atoms { get; set; }
    }

    /// <summary>
    ///      Represents an electronic structure problem.
    /// </summary>
    /// <remarks>
    ///      This type is agnostic with respect to serialization formats, and is a proper subset of
    ///      all currently supported formats such that all current serialization formats should correctly
    ///      serialize and deserialize to this type.
    ///
    ///      Fields and properties that can be missing are explicitly marked as nullable.
    /// </remarks>
    public class ElectronicStructureProblem : IHasMetadata
    {
        #region Metadata
        
        /// <summary>
        ///     Represents any additional metadata about this electronic structure problem.
        /// </summary>
        public Dictionary<string, object> Metadata { get; set; } = new Dictionary<string, object>();

        public Quantity<double> EnergyOffset { get; set; }
        public Quantity<double> CoulombRepulsion { get; set; }

        /// <summary>
        ///     The self-consistent field energy for this electronic structure problem.
        ///     If no energy SCF is provided or known, may be <c>null</c>.
        /// </summary>
        public Quantity<double>? ScfEnergy { get; set; }
        public Quantity<double>? ScfEnergyOffset { get; set; }

        /// <summary>
        ///     The full configuration interaction energy for this electronic structure problem.
        ///     If the FCI energy is known only to be bounded with in an interval,
        ///     the <see cref="BoundedQuantity{TValue}.Value" /> property may be <c>null</c>,
        ///     while if not even a bound is known, this property itself may be <c>null</c>.
        /// </summary>
        public BoundedQuantity<double>? FciEnergy { get; set; }

        public BasisSet? BasisSet { get; set; }

        public int NOrbitals { get; set; }
        public int NElectrons { get; set; }

        public Geometry? Geometry { get; set; }

        #endregion

        
        /// <summary>
        ///     Hamiltonian represented by orbital integrals.
        /// </summary>
        public OrbitalIntegralHamiltonian OrbitalIntegralHamiltonian { get; set; }

        /// <summary>
        ///     Optional collection of trial wavefunctions.
        /// </summary>
        public Dictionary<string, FermionWavefunction<SpinOrbital>>? InitialStates { get; set; }

    }

}

