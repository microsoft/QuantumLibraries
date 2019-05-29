// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.
using System;
using System.Linq;
using System.Collections.Generic;

using Microsoft.Quantum.Simulation.Core;

using Microsoft.Quantum.Chemistry.Broombridge;
using Microsoft.Quantum.Chemistry.OrbitalIntegrals;
using Microsoft.Quantum.Chemistry.Fermion;
using Microsoft.Quantum.Chemistry.Pauli;
using Microsoft.Quantum.Chemistry.QSharpFormat;
using Microsoft.Quantum.Chemistry.JordanWigner;
using Microsoft.Quantum.Chemistry.Generic;


namespace Microsoft.Quantum.Chemistry
{
    /// <summary>
    /// This class contains convenience functions for simulating electronic structure problems.
    /// </summary>
    public static class Convenience
    {
        public class ProblemContainer
        {
            public class Config { }
            //public readonly Config.IndexConvention.LadderType IndexConvention;

            // For now, this only has a few options.
            // Hamiltonian data structures.
            public OrbitalIntegralHamiltonian orbitalIntegralHamiltonian = new OrbitalIntegralHamiltonian();
            public FermionHamiltonian fermionHamiltonian = new FermionHamiltonian();
            public PauliHamiltonian pauliHamiltonian = new PauliHamiltonian();

            // Wavefunction data structure.
            //public Dictionary<string, InputState> wavefunctions = new Dictionary<string, InputState>();

            // QSharp data structure.
            public JordanWignerEncodingData qSharpData;
            
            // Additional data
            public int NOrbitals = 0;
            public int NElectrons = 0;
            public string MiscellaneousInformation;

            /// <summary>
            ///     A label for this particular Hamiltonian.
            ///     Can be used to identify the Hamiltonian out of set
            ///     loaded from the same file.
            /// </summary>
            public string Name { get; set; } = "<unknown>";

        }

        /// <summary>
        ///      Loads a Hamiltonian from integral data represented
        ///      in Broombridge format.
        ///      Please see the <a href="https://docs.microsoft.com/quantum/libraries/chemistry/schema/spec">
        ///      for further details about the
        ///      format parsed by this method.
        /// </summary>
        /// <param name="filename">The name of the file to be loaded.</param>
        /// <returns>
        ///      An instance of <see cref="FermionHamiltonian"/> representing the
        ///      data contained in <paramref name="filename"/>.
        /// </returns>
        public static IEnumerable<FermionHamiltonian> LoadFromBroombridge(
            string filename,
            IndexConvention indexConvention)
        {
            var broombridgeData = Deserializers.DeserializeBroombridge(filename);
            
            // Create electronic structure Hamiltonian
            var fermionHamiltonians = broombridgeData.ProblemDescriptions
                .Select(o => o.OrbitalIntegralHamiltonian
                .ToFermionHamiltonian(indexConvention));

            return fermionHamiltonians;
        }

    }

    public static partial class Extensions
    {
        /// <summary>
        /// Converts an electronic structure problem description
        /// into a format consumable by Q# using default settings.
        /// </summary>
        /// <param name="problem">Input electronic structure problem description.</param>
        /// <param name="state">Selected wavefunction ansatz. This uses the Hartree-Fock state by default.</param>
        /// <param name="indexConvention">Convention for mapping spin-orbit indices to integer indices.</param>
        /// <param name="qubitEncoding">Scheme for mapping fermions to qubits.</param>
        /// <returns></returns>
        public static JordanWignerEncodingData ToQSharpFormat(
            this ProblemDescription problem,
            string state = "",
            IndexConvention indexConvention = IndexConvention.UpDown,
            QubitEncoding qubitEncoding = QubitEncoding.JordanWigner
            )
        {
            var fermionHamiltonian = problem
                .OrbitalIntegralHamiltonian
                .ToFermionHamiltonian(indexConvention);

            var wavefunction = problem.Wavefunctions.ContainsKey(state) ?
                problem.Wavefunctions[state].ToIndexing(indexConvention) :
                fermionHamiltonian.CreateHartreeFockState(problem.NElectrons);

            var pauliHamiltonian = fermionHamiltonian.ToPauliHamiltonian(qubitEncoding);

            var pauliHamiltonianQSharpFormat = pauliHamiltonian.ToQSharpFormat();
            var wavefunctionQSharpFormat = wavefunction.ToQSharpFormat();

            return QSharpFormat.Convert.ToQSharpFormat(pauliHamiltonianQSharpFormat, wavefunctionQSharpFormat);
        }
    }
}