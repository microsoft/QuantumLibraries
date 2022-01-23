using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

using Microsoft.Jupyter.Core;
using Microsoft.Quantum.Chemistry.Broombridge;
using Microsoft.Quantum.Chemistry.Fermion;
using Microsoft.Quantum.Chemistry.LadderOperators;
using Microsoft.Quantum.Chemistry.OrbitalIntegrals;

using Newtonsoft.Json;

namespace Microsoft.Quantum.Chemistry.Magic
{
    /// <summary>
    /// Loads Broombridge electronic structure problem and returns trial wavefunction description.
    /// </summary>
    public class WavefunctionMagic : MagicSymbol
    {
        public WavefunctionMagic()
        {
            this.Name = $"%chemistry.inputstate.load";
            this.Documentation = new Microsoft.Jupyter.Core.Documentation() { Summary = "Loads Broombridge electronic structure problem and returns selected input state." };
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
            [JsonProperty(PropertyName = "file_name")]
            public string FileName { get; set; }

            /// <summary>
            /// A Broombridge ProblemDescription to load the FermionHamiltonian from.
            /// </summary>
            [JsonProperty(PropertyName = "problem_description")]
            public V0_2.ProblemDescription ProblemDescription { get; set; }

            /// <summary>
            /// The IndexConvention to use to generate the Hamiltonian from the ProblemDescription.
            /// </summary>
            [JsonProperty(PropertyName = "index_convention")]
            public IndexConvention IndexConvention { get; set; } = IndexConvention.UpDown;

            /// <summary>
            /// The label of the wavefunction within the ProblemDescription to use.
            /// If no label specified, it will return the Hartree-Fock state.
            /// </summary>
            [JsonProperty(PropertyName = "wavefunction_label")]
            public string WavefunctionLabel { get; set; }
        }

        /// <summary>
        /// Loads a InputState (WaveFunction) from the given Broombridge's ProblemDescription.
        /// The ProblemDescription can come from a .yaml broombridge file or can be passed down as parameter.
        /// If the wavefunctionLabel is specified, it will return the corresponding inputstate from the
        /// ProblemDescription; if the wavefunctionLabel is not specified, then it returns
        /// the Hartree--Fock state.
        /// </summary>
        public Task<ExecutionResult> Run(string input, IChannel channel)
        {
            if (string.IsNullOrWhiteSpace(input))
            {
                channel.Stderr("Please provide the name of a Broombridge file or a problem description to load the fermion Hamiltonian from.");
                return Task.FromResult(ExecuteStatus.Error.ToExecutionResult());
            }

            // Identify the ProblemDescription with the hamiltonian from the arguments.
            var args = JsonConvert.DeserializeObject<Arguments>(input);
            var problemData = SelectProblemDescription(args);

            // Based on the argument, return the Hartree--Fock state or the wavefunction with the given label.
            var wavefunction = (string.IsNullOrEmpty(args.WavefunctionLabel))
                ? problemData.OrbitalIntegralHamiltonian.ToFermionHamiltonian(args.IndexConvention).CreateHartreeFockState(problemData.NElectrons)
                : problemData.Wavefunctions[args.WavefunctionLabel].ToIndexing(args.IndexConvention);

            return Task.FromResult(wavefunction.ToExecutionResult());
        }

        /// <summary>
        /// Selects the ProblemDescription from the given arguments.
        /// If the fileName is specified, it will try to load the Broombridge data from the file
        /// and will use the first ProblemDescription.
        /// </summary>
        protected virtual ProblemDescription SelectProblemDescription(Arguments args)
        {
            if (string.IsNullOrWhiteSpace(args.FileName))
            {
                return Broombridge.ProblemDescription.ProcessRawProblemDescription(args.ProblemDescription);
            }

            // A single file can contain multiple problem descriptions. Let us pick the first one.
            Data broombridge = Deserializers.DeserializeBroombridge(args.FileName);
            return broombridge.ProblemDescriptions.First();
        }
    }
}