// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.QAOA.Jupyter
{
    using System.Threading.Tasks;
    using Microsoft.Jupyter.Core;
    using Microsoft.Quantum.QAOA.QaoaHybrid;
    using Newtonsoft.Json;

    public class HybridQaoaProblemInstanceMagic : MagicSymbol
    {
        public HybridQaoaProblemInstanceMagic()
        {
            this.Name = $"%qaoa.hybridqaoa.create.problem.instance";
            this.Documentation = new Documentation() { Summary = "Prepares a problem instance object that serves as one of arguments to %qaoa.hybridqaoa.run." };
            this.Kind = SymbolKind.Magic;
            this.Execute = this.Run;
        }

        /// <summary>
        /// List of arguments.
        /// </summary>
        public class Arguments
        {
            /// <summary>
            /// Coefficents for one-local Hamiltonian terms.
            /// </summary>
            [JsonProperty(PropertyName = "one_local_hamiltonian_coefficients")]
            public double[] OneLocalHamiltonianCoefficients { get; set; }

            /// <summary>
            /// Coefficents for one-local Hamiltonian terms.
            /// </summary>
            [JsonProperty(PropertyName = "two_local_hamiltonian_coefficients")]
            public double[] TwoLocalHamiltonianCoefficients { get; set; }

            /// <summary>
            /// Size of the combinatorial problem in bits.
            /// </summary>
            [JsonProperty(PropertyName = "problem_size_in_bits")]
            public int ProblemSizeInBits { get; set; }
        }

        /// <summary>
        /// Prepares a ProblemInstance object for a hybrid QAOA.
        /// </summary>
        public async Task<ExecutionResult> Run(string input, IChannel channel)
        {
            if (string.IsNullOrWhiteSpace(input))
            {
                channel.Stderr("Please provide correct arguments");
                return ExecuteStatus.Error.ToExecutionResult();
            }

            var args = JsonConvert.DeserializeObject<Arguments>(input);

            var problemInstance = new ProblemInstance(args.OneLocalHamiltonianCoefficients, args.TwoLocalHamiltonianCoefficients);

            return problemInstance.ToExecutionResult();
        }
    }
}
