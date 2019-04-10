using System.Collections.Generic;
using System.Linq;

using Microsoft.Jupyter.Core;
using Microsoft.Quantum.Chemistry.Broombridge;
using Microsoft.Quantum.Chemistry.Fermion;
using Microsoft.Quantum.Chemistry.LadderOperators;
using Microsoft.Quantum.Chemistry.OrbitalIntegrals;

using Newtonsoft.Json;

namespace Microsoft.Quantum.Chemistry.Magic
{
    /// <summary>
    /// Loads Broombridge electronic structure problem and returns fermion Hamiltonian.
    /// </summary>
    public class WavefunctionMagic : MagicSymbol
    {
        public WavefunctionMagic()
        {
            this.Name = $"%inputstate-load";
            this.Documentation = new Documentation() { Summary = "Loads Broombridge electronic structure problem and returns selected input state." };
            this.Kind = SymbolKind.Magic;
            this.Execute = this.Run;
        }

        /// <summary>
        /// List of arguments
        /// </summary>
        public class Arguments
        {
            /// <summary>
            /// The name of a .yaml file with a broombridge schema.
            /// </summary>
            public string fileName { get; set; }

            /// <summary>
            /// A Broombridge ProblemDescription to load the FermionHamiltonian from.
            /// </summary>
            public Data.ProblemDescription problemDescription { get; set; }

            /// <summary>
            /// The IndexConvention to use to generate the Hamiltonian from the ProblemDescription.
            /// </summary>
            public IndexConvention indexConvention { get; set; } = IndexConvention.UpDown;

            /// <summary>
            /// The label of the wavefunctio within the ProblemDescription to use. 
            /// If no label specified, it will return the Hartree-Fock state.
            /// </summary>
            public string wavefunctionLabel { get; set; }
        }

        /// <summary>
        /// Loads a InputState (WaveFunction) from the given Broombridge's ProblemDescription.
        /// The ProblemDescription can come from a .yaml broombridge file or can be passed down as parameter.
        /// If the wavefunctionLabel is specified, it will return the corresponding inputstate from the
        /// ProblemDescription; if the wavefunctionLabel is not specified, then it returns 
        /// the Hartree--Fock state.
        /// </summary>
        public ExecutionResult Run(string input, IChannel channel)
        {
            if (string.IsNullOrWhiteSpace(input))
            {
                channel.Stderr("Please provide the name of a Broombridge file or a problem description to load the fermion Hamiltonian from.");
                return ExecuteStatus.Error.ToExecutionResult();
            }

            // Identify the ProblemDescription with the hamiltonian from the arguments.
            var args = JsonConvert.DeserializeObject<Arguments>(input);
            var problemData = SelectProblemDescription(args);

            // Based on the argument, return the Hartree--Fock state or the wavefunction with the given label.
            var inputState = (string.IsNullOrEmpty(args.wavefunctionLabel))
                ? problemData.OrbitalIntegralHamiltonian.ToFermionHamiltonian(args.indexConvention).CreateHartreeFockState(problemData.NElectrons)
                : problemData.Wavefunctions[args.wavefunctionLabel].ToIndexing(args.indexConvention);

            return inputState.ToExecutionResult();
        }
        
        /// <summary>
        /// Selects the ProblemDescription from the given arguments.
        /// If the fileName is specified, it will try to load the Broombridge data from the file
        /// and will use the first ProblemDescription, otherwise, the problemDescription in the arguments is used.
        /// </summary>
        protected virtual Data.ProblemDescription SelectProblemDescription(Arguments args)
        {
            if (string.IsNullOrWhiteSpace(args.fileName))
            {
                return args.problemDescription;
            }

            // A single file can contain multiple problem descriptions. Let us pick the first one.
            Data broombridge = Deserializers.DeserializeBroombridge(args.fileName);
            return broombridge.ProblemDescriptions.First();
        }
    }
}