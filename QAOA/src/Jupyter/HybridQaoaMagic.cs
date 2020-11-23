// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Qaoa.Jupyter
{
    using System;
    using System.Threading.Tasks;
    using Microsoft.Jupyter.Core;
    using Microsoft.Quantum.Qaoa.QaoaHybrid;
    using Newtonsoft.Json;

    /// <summary>
    /// This class makes it possible to run the hybrid QAOA with user-defined input parameters in languages other than C# that are supported by this mechanism.
    /// </summary>
    public class HybridQaoaRunMagic : MagicSymbol
    {

        public HybridQaoaRunMagic()
        {
            this.Name = $"%qaoa.hybridqaoa.run";
            this.Documentation = new Documentation()
            {
                Summary = "Runs a hybrid QAOA algorithm with a classical optimizer that chooses input angles. QAOA parameters are provided as JSON input. Initial beta and gamma coefficients are provided by a user."
            };
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
            /// A parameter related to the depth of a QAOA circuit.
            /// </summary>
            [JsonProperty(PropertyName = "NHamiltonianApplications")]
            public int NHamiltonianApplications { get; set; }

            /// <summary>
            /// Description of a combinatorial problem to be solved.
            /// </summary>
            [JsonProperty(PropertyName = "problem_instance")]
            public QaoaProblemInstance QaoaProblemInstance { get; set; }

            /// <summary>
            /// Initial QAOA parameters (beta and gamma angles).
            /// </summary>
            [JsonProperty(PropertyName = "initial_qaoa_parameters")]
            public QaoaParameters InitialQaoaParameters { get; set; }

            /// <summary>
            /// Flag whether optimization should be logged into a file.
            /// </summary>
            [JsonProperty(PropertyName = "should_log")]
            public Boolean ShouldLog { get; set; } = false;
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

            var hybridQaoa = new HybridQaoa(args.NumberOfIterations, args.NHamiltonianApplications, args.QaoaProblemInstance, args.ShouldLog);

            return hybridQaoa.RunOptimization(args.InitialQaoaParameters).ToExecutionResult();
        }
    }

    /// <summary>
    /// This class makes it possible to run the hybrid QAOA with randomly initialized input parameters in languages other than C# that are supported by this mechanism.
    /// </summary>
    public class HybridQaoaWithRandomParametersRunMagic : MagicSymbol
    {

        public HybridQaoaWithRandomParametersRunMagic()
        {
            this.Name = $"%qaoa.hybridqaoa.random.params.run";
            this.Documentation = new Documentation()
            {
                Summary = "Runs a hybrid QAOA algorithm with a classical optimizer that chooses input angles. QAOA parameters are provided as a JSON. Initial beta and gamma parameters are chosen randomly."
            };
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
            /// A parameter related to the depth of a QAOA circuit.
            /// </summary>
            [JsonProperty(PropertyName = "NHamiltonianApplications")]
            public int NHamiltonianApplications { get; set; }

            /// <summary>
            /// Description of a combinatorial problem to be solved.
            /// </summary>
            [JsonProperty(PropertyName = "problem_instance")]
            public QaoaProblemInstance QaoaProblemInstance { get; set; }

            /// <summary>
            /// Number of random starting points in the angles parameters spaces.
            /// </summary>
            [JsonProperty(PropertyName = "number_of_random_starting_points")]
            public int NumberOfRandomStartingPoints { get; set; } = 1;

            /// <summary>
            /// Flag whether optimization should be logged into a file.
            /// </summary>
            [JsonProperty(PropertyName = "should_log")]
            public Boolean ShouldLog { get; set; } = false;
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

            var hybridQaoa = new HybridQaoa(args.NumberOfIterations, args.NHamiltonianApplications, args.QaoaProblemInstance, args.ShouldLog);

            return hybridQaoa.RunOptimization(args.NumberOfRandomStartingPoints).ToExecutionResult();
        }
    }
}
