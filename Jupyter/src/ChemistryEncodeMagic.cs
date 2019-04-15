
using Microsoft.Jupyter.Core;
using Microsoft.Quantum.Chemistry.Fermion;
using Microsoft.Quantum.Chemistry.LadderOperators;
using Microsoft.Quantum.Chemistry.Pauli;
using Microsoft.Quantum.Chemistry.QSharpFormat;

using Newtonsoft.Json;

namespace Microsoft.Quantum.Chemistry.Magic
{
    /// <summary>
    /// Encodes a fermion Hamiltonian and wavefunction ansatz into a format consumable by Q#.
    /// </summary>
    public class ChemistryEncodeMagic : MagicSymbol
    {
        public ChemistryEncodeMagic()
        {
            this.Name = $"%chemistry.encode";
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
            public FermionWavefunction<int> wavefunction { get; set; }
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
            var qSharpWavefunction = args.wavefunction.ToQSharpFormat();
            var qSharpData = Microsoft.Quantum.Chemistry.QSharpFormat.Convert.ToQSharpFormat(qSharpHamiltonian, qSharpWavefunction);

            return qSharpData.ToExecutionResult();
        }
    }
}
