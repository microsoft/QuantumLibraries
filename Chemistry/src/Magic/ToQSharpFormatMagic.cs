using System;
using System.Collections.Generic;
using System.Linq;
using Microsoft.Jupyter.Core;
using Microsoft.Quantum.Chemistry;
using Microsoft.Quantum.IQSharp;
using Newtonsoft.Json;
using Newtonsoft.Json.Linq;

using Microsoft.Quantum.Chemistry.Broombridge;
using Microsoft.Quantum.Chemistry.OrbitalIntegrals;
using Microsoft.Quantum.Chemistry.Fermion;
using Microsoft.Quantum.Chemistry.Pauli;
using Microsoft.Quantum.Chemistry.QSharpFormat;
using Microsoft.Quantum.Chemistry.Generic;
using Microsoft.Quantum.Chemistry;

namespace Magic
{
 

    /// <summary>
    /// "Converts a fermion Hamiltonian and wavefunction ansatz to a format consumable by Q#."
    /// </summary>
    public class ToQSharpFormat : MagicSymbol
    {
        public ToQSharpFormat()
        {
            this.Name = $"%toQSharpFormat";
            this.Documentation = new Documentation() { Summary = "Converts a fermion Hamiltonian to a format consumable by Q#." };
            this.Kind = SymbolKind.Magic;
            this.Execute = this.Run;
        }

        public class Arguments
        {
            public (int[], Dictionary<TermType.Fermion, (int[], double)[]>) serializedHamiltonian { get; set; }
            public (Microsoft.Quantum.Chemistry.StateType, (double, (RaisingLowering, int)[])[]) serializedInputState { get; set; }
        }

        public ExecutionResult Run(string input, IChannel channel)
        {
            var args = Newtonsoft.Json.JsonConvert.DeserializeObject<Arguments>(input);

            var fermionHamiltonian = new FermionHamiltonian(args.serializedHamiltonian);
            var inputState = new InputState(args.serializedInputState);

            // We target a qubit quantum computer, which requires a Pauli representation of the fermion Hamiltonian.
            // A number of mappings from fermions to qubits are possible. Let us choose the Jordan-Wigner encoding.
            PauliHamiltonian pauliHamiltonian = fermionHamiltonian.ToPauliHamiltonian(QubitEncoding.JordanWigner);

            // We now convert this Hamiltonian and a selected state to a format that than be passed onto the QSharp component
            // of the library that implements quantum simulation algorithms.
            var qSharpHamiltonian = pauliHamiltonian.ToQSharpFormat();
            var qSharpWavefunction = inputState.ToQSharpFormat();
            var qSharpData = Microsoft.Quantum.Chemistry.QSharpFormat.Convert.ToQSharpFormat(qSharpHamiltonian, qSharpWavefunction);

            return qSharpData.ToExecutionResult();
        }
    }
}