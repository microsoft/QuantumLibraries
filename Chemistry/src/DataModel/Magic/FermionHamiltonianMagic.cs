using System.Collections.Generic;
using System.Linq;

using Microsoft.Jupyter.Core;
using Microsoft.Quantum.Chemistry.Broombridge;
using Microsoft.Quantum.Chemistry.Fermion;
using Microsoft.Quantum.Chemistry.Generic;
using Microsoft.Quantum.Chemistry.OrbitalIntegrals;
using Microsoft.Quantum.Chemistry.Pauli;
using Microsoft.Quantum.Chemistry.QSharpFormat;

using Newtonsoft.Json;

namespace Microsoft.Quantum.Chemistry.Magic
{
    /// <summary>
    /// Loads Broombridge electronic structure problem and returns fermion Hamiltonian.
    /// </summary>
    public class FermionHamiltonianLoadMagic : MagicSymbol
    {
        public FermionHamiltonianLoadMagic()
        {
            this.Name = $"%fh_load";
            this.Documentation = new Documentation() { Summary = "Loads Broombridge electronic structure problem and returns fermion Hamiltonian." };
            this.Kind = SymbolKind.Magic;
            this.Execute = this.Run;
        }

        public ExecutionResult Run(string input, IChannel channel)
        {
            if (string.IsNullOrWhiteSpace(input))
            {
                channel.Stderr("Please provide the name of a broombridge file to load\n");
                return ExecuteStatus.Error.ToExecutionResult();
            }

            // Deserialize Broombridge from file.
            CurrentVersion.Data broombridge = Deserializers.DeserializeBroombridge(input);

            // A single file can contain multiple problem descriptions. Let us pick the first one.
            CurrentVersion.ProblemDescription problemData = broombridge.ProblemDescriptions.First();
            
            // Electronic structure Hamiltonians are usually represented compactly by orbital integrals. Let us construct
            // such a Hamiltonian from broombridge.
            OrbitalIntegralHamiltonian orbitalIntegralHamiltonian = problemData.ToOrbitalIntegralHamiltonian();

            // We can obtain the full fermion Hamiltonian from the more compact orbital integral representation.
            // This transformation requires us to pick a convention for converting a spin-orbital index to a single integer.
            // Let us pick one according to the formula `integer = 2 * orbitalIndex + spinIndex`.
            FermionHamiltonian fermionHamiltonian = orbitalIntegralHamiltonian.ToFermionHamiltonian(SpinOrbital.IndexConvention.UpDown);
            
            return fermionHamiltonian.ToExecutionResult();
        }
    }

    /// <summary>
    /// Creates an empty fermion Hamiltonian instance.
    /// </summary>
    public class FermionHamiltonianCreateMagic : MagicSymbol
    {
        public FermionHamiltonianCreateMagic()
        {
            this.Name = $"%fh_create";
            this.Documentation = new Documentation() { Summary = "Creates an empty fermion Hamiltonian instance." };
            this.Kind = SymbolKind.Magic;
            this.Execute = this.Run;
        }

        public ExecutionResult Run(string input, IChannel channel)
        {
            // Create empty fermion Hamiltonian instance.
            FermionHamiltonian fermionHamiltonian = new FermionHamiltonian();

            return fermionHamiltonian.ToExecutionResult();
        }
    }

    /// <summary>
    /// Adds terms to a fermion Hamiltonian.
    /// </summary>
    public class FermionHamiltonianAddTermsMagic : MagicSymbol
    {
        public FermionHamiltonianAddTermsMagic()
        {
            this.Name = $"%fh_add_terms";
            this.Documentation = new Documentation() { Summary = "Adds terms to a fermion Hamiltonian." };
            this.Kind = SymbolKind.Magic;
            this.Execute = this.Run;
        }

        /// <summary>
        /// The list of arguments
        /// </summary>
        public class Arguments
        {
            /// <summary>
            /// The fermion hamiltonian to add terms to.
            /// </summary>
            public FermionHamiltonian hamiltonian { get; set; }

            /// <summary>
            /// The list of terms and their coefficient to add.
            /// </summary>
            public List<(HermitianFermionTerm, double)> fermionTerms { get; set; }
        }

        public ExecutionResult Run(string input, IChannel channel)
        {
            var args = JsonConvert.DeserializeObject<Arguments>(input);
            args.hamiltonian.AddRange(args.fermionTerms.Select(t => (t.Item1, t.Item2.ToDoubleCoeff())));
            
            return args.hamiltonian.ToExecutionResult();
        }
    }

    /// <summary>
    /// Encodes a fermion Hamiltonian and wavefunction ansatz into a format consumable by Q#.
    /// </summary>
    public class FermionHamiltonianEncodeMagic : MagicSymbol
    {
        public FermionHamiltonianEncodeMagic()
        {
            this.Name = $"%fh_encode";
            this.Documentation = new Documentation() { Summary = "Encodes a fermion Hamiltonian to a format consumable by Q#." };
            this.Kind = SymbolKind.Magic;
            this.Execute = this.Run;
        }

        public class Arguments
        {
            /// <summary>
            /// The fermion hamiltonian.
            /// </summary>
            public FermionHamiltonian hamiltonian { get; set; }

            /// <summary>
            /// The input state.
            /// </summary>
            public InputState inputState { get; set; }
        }

        public ExecutionResult Run(string input, IChannel channel)
        {
            var args = JsonConvert.DeserializeObject<Arguments>(input);

            // We target a qubit quantum computer, which requires a Pauli representation of the fermion Hamiltonian.
            // A number of mappings from fermions to qubits are possible. Let us choose the Jordan-Wigner encoding.
            PauliHamiltonian pauliHamiltonian = args.hamiltonian.ToPauliHamiltonian(QubitEncoding.JordanWigner);

            // We now convert this Hamiltonian and a selected state to a format that than be passed onto the QSharp component
            // of the library that implements quantum simulation algorithms.
            var qSharpHamiltonian = pauliHamiltonian.ToQSharpFormat();
            var qSharpWavefunction = args.inputState.ToQSharpFormat();
            var qSharpData = Microsoft.Quantum.Chemistry.QSharpFormat.Convert.ToQSharpFormat(qSharpHamiltonian, qSharpWavefunction);

            return qSharpData.ToExecutionResult();
        }
    }
}