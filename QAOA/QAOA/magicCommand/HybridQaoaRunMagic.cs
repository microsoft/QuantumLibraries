using System;
using System.Threading.Tasks;
using Microsoft.Jupyter.Core;
using Newtonsoft.Json;
using Quantum.QAOA;

namespace QAOA.magicCommand
{
    public class HybridQaoaRunMagic : MagicSymbol
    {

        public HybridQaoaRunMagic()
        {
            this.Name = $"%standard.hybridqaoa.run";
            this.Documentation = new Documentation() { Summary = "Runs a hybrid QAOA algorithm with a classical optimizer that chooses input angles. QAOA parameters are provided as a json" };
            this.Kind = SymbolKind.Magic;
            this.Execute = this.Run;
        }

        /// <summary>
        /// List of arguments
        /// </summary>
        public class Arguments
        {
            /// <summary>
            /// Number of iterations in the fidelity sampling.
            /// </summary>
            [JsonProperty(PropertyName = "number_of_iterations")]
            public int NumberOfIterations { get; set; }

            /// <summary>
            /// Depth of a QAOA circuit.
            /// </summary>
            [JsonProperty(PropertyName = "p")]
            public int p { get; set; }

            /// <summary>
            /// Description of a combinatorial problem to be solved.
            /// </summary>
            [JsonProperty(PropertyName = "problem_instance")]
            public ProblemInstance ProblemInstance { get; set; }

            /// <summary>
            /// Number of random starting points in the angles parameters spaces.
            /// </summary>
            [JsonProperty(PropertyName = "number_of_random_starting_points")]
            public int NumberOfRandomStartingPoints { get; set; } = 1;

            /// <summary>
            /// Initial beta angles.
            /// </summary>
            [JsonProperty(PropertyName = "initial_beta")]
            public Double[] InitialBeta { get; set; } = null;

            /// <summary>
            /// Initial gamma angles.
            /// </summary>
            [JsonProperty(PropertyName = "initial_gamma")]
            public Double[] InitialGamma { get; set; } = null;
        }

        /// <summary>
        /// Runs a hybrid QAOA.
        /// </summary>
        public async Task<ExecutionResult> Run(string input, IChannel channel)
        {
            if (string.IsNullOrWhiteSpace(input))
            {
                channel.Stderr("Please provide correct arguments");
                return ExecuteStatus.Error.ToExecutionResult();
            }

            var args = JsonConvert.DeserializeObject<Arguments>(input);


            HybridQaoa hybridQaoa = new HybridQaoa(args.NumberOfIterations, args.p, args.ProblemInstance); //deal with optional params

            return hybridQaoa.runOptimization().ToExecutionResult();
        }


    }
}
