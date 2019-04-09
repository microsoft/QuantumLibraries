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

namespace Magic
{
    // Expose function that creates fermion Hamiltonian froom  Broombridge
    // Expose function to create empty FermionHamiltonian
    // Expose function to populate FermionHamiltonian

    // Python: 1) Make empty FermionHamiltonian -- this returns serialization of empty Hamiltonian
    //          2) Python AddTerms to the Hamiltonian -- this passes the Ham back to c#, and returns one populated with terms. Note: Make SetTerms for multiple terms.
    //          3) Python .ToQSharpFormat -- back to C# and back out.

    // Do the same for quantum states.


    /// <summary>
    /// Loads Broombridge electronic structure problem and returns fermion Hamiltonian.
    /// </summary>
    public class FermionHamiltonianFromBroombridge : MagicSymbol
    {
        public FermionHamiltonianFromBroombridge()
        {
            this.Name = $"%fermionHamiltonianFromBroombridge";
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


            // TODO: Implement serialization of fermion Hamiltonian first.
            var fermionHamiltonianData = Newtonsoft.Json.JsonConvert.SerializeObject(fermionHamiltonian);
            
            return fermionHamiltonianData.ToExecutionResult();
        }
    }

    /// <summary>
    /// "Creates an empty fermion Hamiltonian instance." 
    /// </summary>
    public class CreateNewFermionHamiltonian : MagicSymbol
    {
        public CreateNewFermionHamiltonian()
        {
            this.Name = $"%createNewFermionHamiltonian";
            this.Documentation = new Documentation() { Summary = "Creates an empty fermion Hamiltonian instance." };
            this.Kind = SymbolKind.Magic;
            this.Execute = this.Run;
        }

        public ExecutionResult Run(string input, IChannel channel)
        {
            // Create empty fermion Hamiltonian instance.
            FermionHamiltonian fermionHamiltonian = new FermionHamiltonian();

            // TODO: Implement serialization of fermion Hamiltonian first.
            var fermionHamiltonianData = Newtonsoft.Json.JsonConvert.SerializeObject(fermionHamiltonian);

            return fermionHamiltonianData.ToExecutionResult();
        }
    }

    /// <summary>
    /// "Adds terms to a fermion Hamiltonian."
    /// </summary>
    public class AddTermsToFermionHamiltonian : MagicSymbol
    {
        public AddTermsToFermionHamiltonian()
        {
            this.Name = $"%addTermsToFermionHamiltonian";
            this.Documentation = new Documentation() { Summary = "Adds terms to a fermion Hamiltonian." };
            this.Kind = SymbolKind.Magic;
            this.Execute = this.Run;
        }

        public class Arguments
        {
            public Arguments()
            {
                this.fermionHamiltonian = new FermionHamiltonian();
                this.fermionTerms = new List<(HermitianFermionTerm, double)>();
            }

            public FermionHamiltonian fermionHamiltonian { get; set; }
            public List<(HermitianFermionTerm, double)> fermionTerms { get; set; }
        }

        public ExecutionResult Run(string input, IChannel channel)
        {
            var args = Newtonsoft.Json.JsonConvert.DeserializeObject<Arguments>(input);

            args.fermionHamiltonian.AddRange(args.fermionTerms.Select(o => (o.Item1, o.Item2.ToDoubleCoeff())));

            // TODO: Implement serialization of fermion Hamiltonian first.
            var fermionHamiltonianData = Newtonsoft.Json.JsonConvert.SerializeObject(args.fermionHamiltonian);

            return fermionHamiltonianData.ToExecutionResult();
        }
    }

    /// <summary>
    /// "Converts a fermion Hamiltonian and wavefunction ansatz to a format consumable by Q#."
    /// </summary>
    public class ToQSharpFormatFromFermionHamiltonian : MagicSymbol
    {
        public ToQSharpFormatFromFermionHamiltonian()
        {
            this.Name = $"%toQSharpFormatFromFermionHamiltonian";
            this.Documentation = new Documentation() { Summary = "Converts a fermion Hamiltonian and wavefunction ansatz to a format consumable by Q#." };
            this.Kind = SymbolKind.Magic;
            this.Execute = this.Run;
        }

        public class Arguments
        {
            public Arguments()
            {
                this.fermionHamiltonian = new FermionHamiltonian();
                this.inputState = new InputState();
            }

            public FermionHamiltonian fermionHamiltonian { get; set; }
            public InputState inputState { get; set; }
        }

        public ExecutionResult Run(string input, IChannel channel)
        {
            var args = Newtonsoft.Json.JsonConvert.DeserializeObject<Arguments>(input);

            // We target a qubit quantum computer, which requires a Pauli representation of the fermion Hamiltonian.
            // A number of mappings from fermions to qubits are possible. Let us choose the Jordan-Wigner encoding.
            PauliHamiltonian pauliHamiltonian = args.fermionHamiltonian.ToPauliHamiltonian(QubitEncoding.JordanWigner);

            // We now convert this Hamiltonian and a selected state to a format that than be passed onto the QSharp component
            // of the library that implements quantum simulation algorithms.
            var qSharpHamiltonian = pauliHamiltonian.ToQSharpFormat();
            var qSharpWavefunction = args.inputState.ToQSharpFormat();
            var qSharpData = Microsoft.Quantum.Chemistry.QSharpFormat.Convert.ToQSharpFormat(qSharpHamiltonian, qSharpWavefunction);

            return qSharpData.ToExecutionResult();
        }
    }
}