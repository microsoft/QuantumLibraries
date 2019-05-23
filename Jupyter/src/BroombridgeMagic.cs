using Microsoft.Jupyter.Core;
using Microsoft.Quantum.Chemistry.Broombridge;

namespace Microsoft.Quantum.Chemistry.Magic
{
    /// <summary>
    /// Jupyter Magic that loads and returns 
    /// Broombridge electronic structure problem representation from a given .yaml file.
    /// </summary>
    public class BroombridgeMagic : MagicSymbol
    {
        /// <summary>
        /// Default constructor.
        /// </summary>
        public BroombridgeMagic()
        {
            this.Name = $"%chemistry.broombridge";
            this.Documentation = new Documentation() { Summary = "Loads and returns Broombridge electronic structure problem representation from a given .yaml file." };
            this.Kind = SymbolKind.Magic;
            this.Execute = this.Run;
        }

        /// <summary>
        /// Loads the broombridge data from the given .yaml file and returns it.
        /// </summary>
        public ExecutionResult Run(string input, IChannel channel)
        {
            if (string.IsNullOrWhiteSpace(input))
            {
                channel.Stderr("Please provide the name of a broombridge file to load\n");
                return ExecuteStatus.Error.ToExecutionResult();
            }

            var yamlData = Deserializers.DeserializeBroombridge(input);
            return yamlData.ToExecutionResult();
        }
    }    
}
