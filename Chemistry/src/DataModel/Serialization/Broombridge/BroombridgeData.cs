// Copyright (c) Microsoft Corporation. All rights reserved.
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
                Wavefunctions = FromBroombridgeV0_2(problem.InitialStates)
            };
            return problemDescription;
        }

        internal static Dictionary<string, Fermion.FermionWavefunction<OrbitalIntegrals.SpinOrbital>>
            FromBroombridgeV0_2(
                IEnumerable<V0_2.State>? initialStates
            )
        {
            var wavefunctions = new Dictionary<string, Fermion.FermionWavefunction<OrbitalIntegrals.SpinOrbital>>();
            foreach (var initialState in initialStates ?? new List<V0_2.State>())
            {
                var finalState = new FermionWavefunction<SpinOrbital>();

                var (method, energy, outputState) = V0_2.ToWavefunction(initialState);

                var setMethod = V0_2.ParseInitialStateMethod(initialState.Method);
                var setEnergy = energy;

                if (setMethod == StateType.SparseMultiConfigurational)
                {
                    var mcfData = (SparseMultiCFWavefunction<SpinOrbital>)outputState;

                    finalState = new FermionWavefunction<SpinOrbital>(
                        mcfData.Excitations
                            .Select(o => (
                                o.Key.ToIndices().ToArray(),
                                o.Value.Real
                            )));
                }
                else if (setMethod == StateType.UnitaryCoupledCluster)
                {
                    var uccData = (UnitaryCCWavefunction<SpinOrbital>)outputState;

                    var reference = uccData.Reference;

                    var excitations = uccData.Excitations;

                    finalState = new FermionWavefunction<SpinOrbital>(
                        reference.ToIndices(),
                        excitations
                            .Select(o => (
                                o.Key.ToIndices().ToArray(),
                                o.Value.Real
                            ))
                        );
                }
                else
                {
                    throw new System.ArgumentException($"Wavefunction type {setMethod} not recognized");
                }

                finalState.Method = setMethod;
                finalState.Energy = setEnergy;

                wavefunctions.Add(initialState.Label, finalState);
            }
            return wavefunctions;
        }
    }

}