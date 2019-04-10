using System.Collections.Generic;
using System.Linq;

using Microsoft.Jupyter.Core;
using Microsoft.Quantum.Chemistry.Broombridge;
using Microsoft.Quantum.Chemistry.Fermion;
using Microsoft.Quantum.Chemistry.Generic;
using Microsoft.Quantum.Chemistry.OrbitalIntegrals;

using Newtonsoft.Json;

namespace Microsoft.Quantum.Chemistry.Magic
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
    public class InputStateFromBroombridge : MagicSymbol
    {
        public InputStateFromBroombridge()
        {
            this.Name = $"%input_encode";
            this.Documentation = new Documentation() { Summary = "Loads Broombridge electronic structure problem and returns selected input state." };
            this.Kind = SymbolKind.Magic;
            this.Execute = this.Run;
        }

        public class Arguments
        {
            public Arguments()
            {
                this.wavefunctionLabel = "Greedy";

            }

            public string broombridge { get; set; }

            public string wavefunctionLabel { get; set; }
        }

        public ExecutionResult Run(string input, IChannel channel)
        {
            if (string.IsNullOrWhiteSpace(input))
            {
                channel.Stderr("Please provide the name of a broombridge file to load\n");
                return ExecuteStatus.Error.ToExecutionResult();
            }

            var args = JsonConvert.DeserializeObject<Arguments>(input);

            // Deserialize Broombridge from file.
            CurrentVersion.Data broombridge = Deserializers.DeserializeBroombridge(args.broombridge);

            // A single file can contain multiple problem descriptions. Let us pick the first one.
            CurrentVersion.ProblemDescription problemData = broombridge.ProblemDescriptions.First();

            #region Create wavefunction Ansatzes
            // A list of trial wavefunctions can be provided in the Broombridge file. For instance, the wavefunction
            // may be a single-reference Hartree--Fock state, a multi-reference state, or a unitary coupled-cluster state.
            Dictionary<string, InputState> inputStates = problemData.ToWavefunctions(SpinOrbital.IndexConvention.UpDown);

            InputState inputState = new InputState();

            // If no states are provided, use the Hartree--Fock state.
            if (inputStates.Count() != 0)
            {
                OrbitalIntegralHamiltonian orbitalIntegralHamiltonian = problemData.ToOrbitalIntegralHamiltonian();
                FermionHamiltonian fermionHamiltonian = orbitalIntegralHamiltonian.ToFermionHamiltonian(SpinOrbital.IndexConvention.UpDown);
                inputState = fermionHamiltonian.GreedyStatePreparation(problemData.NElectrons);
            }
            else
            {
                inputState = inputStates[args.wavefunctionLabel];
            }
            #endregion

            return inputState.ToExecutionResult();
        }
    }

    /// <summary>
    /// "Creates an empty fermion Hamiltonian instance." 
    /// </summary>
    public class InputStateEncodeMagic: MagicSymbol
    {
        public InputStateEncodeMagic()
        {
            this.Name = $"%input_encode";
            this.Documentation = new Documentation() { Summary = "Converts an input state to a format consumable by Q#." };
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
}