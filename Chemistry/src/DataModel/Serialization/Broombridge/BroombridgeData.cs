// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

#nullable enable

using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Text;
using System.Linq;
using System.Text.RegularExpressions;
using YamlDotNet.Core;
using YamlDotNet.Serialization;


using Microsoft.Quantum.Chemistry.Fermion;
using Microsoft.Quantum.Chemistry.OrbitalIntegrals;

namespace Microsoft.Quantum.Chemistry.Broombridge
{
    /// <summary>
    /// Latest Broombridge format.
    /// </summary>
    [Obsolete(
        "Please use collections of ElectronicStructureProblem instead.",
        error: false
    )]
    public class Data
    {
        /// <summary>
        /// Constructor for empty Broombridge data structure.
        /// </summary>
        public Data()
        {
            ProblemDescriptions = new List<ProblemDescription>();
        }

        /// <summary>
        /// Raw deserialized Broombridge data.
        /// </summary>
        public V0_2.Data Raw { get; set; }

        // Root of Broombridge data structure

        /// <summary>
        /// URL to schema defining this version of Broombridge.
        /// </summary>
        public string Schema { get; set; }

        /// <summary>
        /// Broombridge instance version number
        /// </summary>
        public VersionNumber VersionNumber { get; set; }

        /// <summary>
        /// Collection of electronic structure problems.
        /// </summary>
        public IEnumerable<ProblemDescription> ProblemDescriptions { get; set; }

        
        /// <summary>
        /// Deserialized Broombridge data
        /// </summary>
        /// <param name="broombridgeV0_2">Broombridge data structure.</param>
        internal Data(Broombridge.V0_2.Data broombridgeV0_2)
        {
            Raw = broombridgeV0_2;
            Schema = broombridgeV0_2.Schema;
            VersionNumber = VersionNumber.v0_2;

            ProblemDescriptions = Raw.ProblemDescriptions.Select(problem => ProblemDescription.ProcessRawProblemDescription(problem));
        }

    }

    /// <summary>
    /// Electronic structure problem instance.
    /// </summary>
    [Obsolete(
        "Please use ElectronicStructureProblem instead.",
        error: false
    )]
    public struct ProblemDescription
    {
        /// <summary>
        /// Identity term of the Hamilonian.
        /// </summary>
        public double EnergyOffset { get; set; }

        /// <summary>
        /// Number of orbitals.
        /// </summary>
        public int NOrbitals { get; set; }

        /// <summary>
        /// Number of electroncs.
        /// </summary>
        public int NElectrons { get; set; }

        /// <summary>
        /// Hamiltonian represented by orbital integrals.
        /// </summary>
        public OrbitalIntegralHamiltonian OrbitalIntegralHamiltonian { get; set; }

        /// <summary>
        /// Collection of trial wavefunctions.
        /// </summary>
        public Dictionary<string, FermionWavefunction<SpinOrbital>> Wavefunctions { get; set; }

        /// <summary>
        /// Converts the raw problem description in V0_2 Broombridge into
        /// in internal data format.
        /// </summary>
        /// <param name="problem">Problem description to be converted</param>
        /// <returns>The internal problem description data structure.</returns>
        public static ProblemDescription ProcessRawProblemDescription(Broombridge.V0_2.ProblemDescription problem)
        {
            var problemDescription = new ProblemDescription
            {
                EnergyOffset = problem.EnergyOffset.Value + problem.CoulombRepulsion.Value,
                NElectrons = problem.NElectrons,
                NOrbitals = problem.NOrbitals,
                OrbitalIntegralHamiltonian = V0_2.ToOrbitalIntegralHamiltonian(problem),
                Wavefunctions = problem.InitialStates?.FromBroombridgeV0_2() ?? new Dictionary<string, FermionWavefunction<SpinOrbital>>()
            };
            return problemDescription;
        }
    }

}