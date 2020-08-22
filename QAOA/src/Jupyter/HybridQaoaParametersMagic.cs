using System.Threading.Tasks;
using Microsoft.Jupyter.Core;
using Microsoft.Quantum.QAOA.QaoaHybrid;
using Newtonsoft.Json;

namespace Microsoft.Quantum.QAOA.Jupyter
{
    public class HybridQaoaParametersMagic : MagicSymbol
    {
        public HybridQaoaParametersMagic()
        {
            this.Name = $"%qaoa.hybridqaoa.create.parameters";
            this.Documentation = new Documentation() { Summary = "Prepares a QAOA parameters object that serves as one of arguments to %qaoa.hybridqaoa.run." };
            this.Kind = SymbolKind.Magic;
            this.Execute = this.Run;
        }

        /// <summary>
        /// List of arguments.
        /// </summary>
        public class Arguments
        {
            /// <summary>
            /// Beta QAOA coefficients.
            /// </summary>
            [JsonProperty(PropertyName = "beta")]
            public double[] Beta { get; set; }

            /// <summary>
            /// Gamma QAOA coefficients.
            /// </summary>
            [JsonProperty(PropertyName = "gamma")]
            public double[] Gamma { get; set; }
        }

        /// <summary>
        /// Prepares a QaoaParameters object for a hybrid QAOA.
        /// </summary>
        public async Task<ExecutionResult> Run(string input, IChannel channel)
        {
            if (string.IsNullOrWhiteSpace(input))
            {
                channel.Stderr("Please provide correct arguments");
                return ExecuteStatus.Error.ToExecutionResult();
            }

            var args = JsonConvert.DeserializeObject<Arguments>(input);

            var qaoaParameters = new QaoaParameters(args.Beta, args.Gamma);

            return qaoaParameters.ToExecutionResult();
        }
    }
}
