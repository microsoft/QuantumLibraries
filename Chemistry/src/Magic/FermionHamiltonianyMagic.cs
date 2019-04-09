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
            var fermionHamiltonianData = Newtonsoft.Json.JsonConvert.SerializeObject(fermionHamiltonian.SerializationFormat());
            
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
            var fermionHamiltonianData = Newtonsoft.Json.JsonConvert.SerializeObject(fermionHamiltonian.SerializationFormat());

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
            public (int[], Dictionary<TermType.Fermion, (int[], double)[]>) serializedHamiltonian { get; set; }
            public List<(int[], double)> fermionTerms { get; set; }
        }

        public ExecutionResult Run(string input, IChannel channel)
        {
            
            var args = Newtonsoft.Json.JsonConvert.DeserializeObject<Arguments>(input);

            var fermionHamiltonian = new FermionHamiltonian(args.serializedHamiltonian);

            fermionHamiltonian.AddRange(args.fermionTerms.Select(o => (new HermitianFermionTerm(o.Item1), o.Item2.ToDoubleCoeff())));
            
            var fermionHamiltonianData = Newtonsoft.Json.JsonConvert.SerializeObject(fermionHamiltonian.SerializationFormat());

            return fermionHamiltonianData.ToExecutionResult();
        }
    }
    
}