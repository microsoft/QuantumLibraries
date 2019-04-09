// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using System.Linq;
using System.Collections.Generic;

using Microsoft.Quantum.Chemistry.Broombridge;
using Microsoft.Quantum.Chemistry.OrbitalIntegrals;
using Microsoft.Quantum.Chemistry.Fermion;
using Microsoft.Quantum.Chemistry.Pauli;
using Microsoft.Quantum.Chemistry.QSharpFormat;
using Microsoft.Quantum.Chemistry.Generic;

namespace Microsoft.Quantum.Chemistry
{
    public class Workflows
    {
        /// <summary>
        /// Sample implementation of end-to-end electronic structure problem simulation. 
        /// </summary>
        /// <param name="filename"></param>
        public static void SampleWorkflow(
            string filename,
            string wavefunctionLabel,
            SpinOrbital.IndexConvention indexConvention
            )
        {
            // Deserialize Broombridge from file.
            CurrentVersion.Data broombridge = Deserializers.DeserializeBroombridge(filename);

            // A single file can contain multiple problem descriptions. Let us pick the first one.
            CurrentVersion.ProblemDescription problemData = broombridge.ProblemDescriptions.First();

            #region Create electronic structure Hamiltonian
            // Electronic structure Hamiltonians are usually represented compactly by orbital integrals. Let us construct
            // such a Hamiltonian from broombridge.
            OrbitalIntegralHamiltonian orbitalIntegralHamiltonian = problemData.ToOrbitalIntegralHamiltonian();

            // We can obtain the full fermion Hamiltonian from the more compact orbital integral representation.
            // This transformation requires us to pick a convention for converting a spin-orbital index to a single integer.
            // Let us pick one according to the formula `integer = 2 * orbitalIndex + spinIndex`.
            FermionHamiltonian fermionHamiltonian = orbitalIntegralHamiltonian.ToFermionHamiltonian(indexConvention);

            // We target a qubit quantum computer, which requires a Pauli representation of the fermion Hamiltonian.
            // A number of mappings from fermions to qubits are possible. Let us choose the Jordan-Wigner encoding.
            PauliHamiltonian pauliHamiltonian = fermionHamiltonian.ToPauliHamiltonian(QubitEncoding.JordanWigner);
            #endregion

            #region Create wavefunction Ansatzes
            // A list of trial wavefunctions can be provided in the Broombridge file. For instance, the wavefunction
            // may be a single-reference Hartree--Fock state, a multi-reference state, or a unitary coupled-cluster state.
            Dictionary<string, InputState> inputStates = problemData.ToWavefunctions(indexConvention);

            // If no states are provided, use the Hartree--Fock state.
            InputState inputState = inputStates.Count() != 0
                ? inputStates[wavefunctionLabel] : fermionHamiltonian.GreedyStatePreparation(problemData.NElectrons);
            #endregion

            #region Pipe to QSharp and simulate
            // We now convert this Hamiltonian and a selected state to a format that than be passed onto the QSharp component
            // of the library that implements quantum simulation algorithms.
            var qSharpHamiltonian = pauliHamiltonian.ToQSharpFormat();
            var qSharpWavefunction = inputStates[wavefunctionLabel].ToQSharpFormat();
            var qSharpData = QSharpFormat.Convert.ToQSharpFormat(qSharpHamiltonian, qSharpWavefunction);
            #endregion
        }
    }

    /// <summary>
    /// Configuration settings for modifying chemistry library behavior.
    /// </summary>
    public class Config
    {

        /// <summary>
        /// Construct default configuration.
        /// </summary>
        /// <returns>Default configuration class.</returns>
        public static Config Default()
        {
            return new Config();
        }

        /// <summary>
        /// Default configuration constructor;
        /// </summary>
        public Config()
        {
            UseIndexConvention = DefaultIndexConvention;
            QubitEncoding UseQubitEncoding = DefaultQubitEncoding;
            UseTruncationThreshold = DefaultTruncationThreshold;
        }

        /// <summary>
        /// Choose indexing convention from spin-orbital index to an integer.
        /// </summary>
        public SpinOrbital.IndexConvention UseIndexConvention;

        /// <summary>
        /// Threshold below which to truncate Hamiltonian coefficients.
        /// </summary>
        public double UseTruncationThreshold;

        /// <summary>
        /// Chose mapping from fermions operators to Pauli operatos.
        /// </summary>
        public QubitEncoding UseQubitEncoding;

        // Default settings
        public const SpinOrbital.IndexConvention DefaultIndexConvention = SpinOrbital.IndexConvention.UpDown;
        public const QubitEncoding DefaultQubitEncoding = QubitEncoding.JordanWigner;
        public const double DefaultTruncationThreshold = 1e-8;
        

    }
        
}
 