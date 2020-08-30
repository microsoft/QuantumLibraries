// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Qaoa.Jupyter
{
    using System.Threading.Tasks;
    using Microsoft.Jupyter.Core;
    using Microsoft.Quantum.Qaoa.QaoaHybrid;
    using Newtonsoft.Json;

    /// <summary>
    /// This class makes it possible to create a problem instance for the hybrid QAOA in languages other than C# that are supported by this mechanism.
    /// </summary>
    public class HybridQaoaProblemInstanceMagic : MagicSymbol
    {
        public HybridQaoaProblemInstanceMagic()
        {
            this.Name = $"%qaoa.hybridqaoa.create.problem.instance";
            this.Documentation = new Documentation()
            {
                Summary = "Prepares a problem instance object that serves as one of arguments to %qaoa.hybridqaoa.run."
            };
            this.Kind = SymbolKind.Magic;
            this.Execute = this.Run;
        }

        /// <summary>
        /// List of arguments.
        /// </summary>
        public class Arguments
        {
            /// <summary>
            /// Coefficients for one-local Hamiltonian terms. Assuming that a solution to a combinatorial optimization problem can be encoded into n bits (which then corresponds to an encoding into n qubits), this array must be of length n. The i-th coefficient in an array corresponds to the coefficient of a term \sigma_i^z.
            /// </summary>
            [JsonProperty(PropertyName = "one_local_hamiltonian_coefficients")]
            public double[] OneLocalHamiltonianCoefficients { get; set; }

            /// <summary>
            /// Coefficients for two-local Hamiltonian terms. Assuming that a solution to a combinatorial optimization problem can be encoded into n bits (which then corresponds to an encoding into n qubits), this array must be of length n^2. The (i*n+j)-th coefficient in an array corresponds to the coefficient of a term \sigma_i^z\sigma_j^z (other coefficients in an array can take any value of type double).
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
        /// Prepares a QaoaProblemInstance object for a hybrid QAOA.
        /// </summary>
        public async Task<ExecutionResult> Run(string input, IChannel channel)
        {
            if (string.IsNullOrWhiteSpace(input))
            {
                channel.Stderr("Please provide correct arguments");
                return ExecuteStatus.Error.ToExecutionResult();
            }

            var args = JsonConvert.DeserializeObject<Arguments>(input);

            var problemInstance = new QaoaProblemInstance(args.OneLocalHamiltonianCoefficients, args.TwoLocalHamiltonianCoefficients);

            return problemInstance.ToExecutionResult();
        }
    }
}
