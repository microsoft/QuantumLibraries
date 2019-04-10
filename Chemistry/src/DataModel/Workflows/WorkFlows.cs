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
            IndexConvention indexConvention
            )
        {
            // Deserialize Broombridge from file.
            Data broombridge = Deserializers.DeserializeBroombridge(filename);

            // A single file can contain multiple problem descriptions. Let us pick the first one.
            Data.ProblemDescription problemData = broombridge.ProblemDescriptions.First();

            #region Create electronic structure Hamiltonian
            // Electronic structure Hamiltonians are usually represented compactly by orbital integrals. Let us construct
            // such a Hamiltonian from broombridge.
            OrbitalIntegralHamiltonian orbitalIntegralHamiltonian = problemData.OrbitalIntegralHamiltonian;

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
            // In this case, Broombridge indexes the fermion operators with spin-orbitals instead of integers. 
            Dictionary<string, FermionWavefunction<SpinOrbital>> inputStates = problemData.Wavefunctions;

            // If no states are provided, use the Hartree--Fock state.
            // As fermion operators the fermion Hamiltonian are already indexed by, we now apply the desired
            // spin-orbital -> integer indexing convention.
            FermionWavefunction<int> inputState = inputStates[wavefunctionLabel].ToIndexing(indexConvention);
            
            //Data.State inputState = inputStates.Count() != 0
            //    ? inputStates[wavefunctionLabel] : fermionHamiltonian.GreedyStatePreparation(problemData.NElectrons);
            #endregion

            #region Pipe to QSharp and simulate
            // We now convert this Hamiltonian and a selected state to a format that than be passed onto the QSharp component
            // of the library that implements quantum simulation algorithms.
            var qSharpHamiltonian = pauliHamiltonian.ToQSharpFormat();
            var qSharpWavefunction = inputState.ToQSharpFormat();
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
            UseIndexConvention = IndexConvention.UpDown;
            QubitEncoding UseQubitEncoding = QubitEncoding.JordanWigner;
            UseTruncationThreshold = 1e-8;
        }

        /// <summary>
        /// Choose indexing convention from spin-orbital index to an integer.
        /// </summary>
        public IndexConvention UseIndexConvention;

        /// <summary>
        /// Threshold below which to truncate Hamiltonian coefficients.
        /// </summary>
        public double UseTruncationThreshold;

        /// <summary>
        /// Chose mapping from fermions operators to Pauli operatos.
        /// </summary>
        public QubitEncoding UseQubitEncoding;
        
        

    }
        
}
 