using System;
using System.IO;
using System.Linq;
using Microsoft.Jupyter.Core;
//using Microsoft.Quantum.Chemistry;
//using Microsoft.Quantum.Simulation.Core;
//using YamlDotNet.Serialization;

namespace Magic
{
    public class BroombridgeMagic //: MagicSymbol
    {

        /*
        public BroombridgeMagic()
        {
            this.Name = $"%broombridge";
            this.Documentation = new Documentation() { Summary = "Loads and returns broombridge hamiltonian representation from a given .yaml file." };
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

            using (var reader = File.OpenText(input))
            {
                var deserializer = new DeserializerBuilder().Build();
                var yamlData = deserializer.Deserialize<IntegralDataSchema>(reader);

                return yamlData.ToExecutionResult();
            }
        }*/
    }
}