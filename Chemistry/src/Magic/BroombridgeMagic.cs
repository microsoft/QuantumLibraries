using System;
using System.IO;
using System.Linq;
using Microsoft.Jupyter.Core;
using Microsoft.Quantum.Chemistry;
using Microsoft.Quantum.Simulation.Core;
using YamlDotNet.Serialization;
using Microsoft.Quantum.Chemistry.Broombridge;

namespace Magic
{
    public class DeserializeBroombridge : MagicSymbol
    {
        public DeserializeBroombridge()
        {
            this.Name = $"%broombridge";
            this.Documentation = new Documentation() { Summary = "Loads and returns Broombridge electronic structure problem representation from a given .yaml file." };
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

            var yamlData = Deserializers.DeserializeBroombridge(input);
            return yamlData.ToExecutionResult();
        }
    }
    
}