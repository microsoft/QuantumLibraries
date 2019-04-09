using System;
using System.Collections.Generic;
using System.Linq;
using Microsoft.Jupyter.Core;
using Microsoft.Quantum.Chemistry;
using Microsoft.Quantum.IQSharp;
using Newtonsoft.Json;
using Newtonsoft.Json.Linq;

namespace Magic
{
    public class JWEncodeMagic : MagicSymbol
    {
        public class Arguments
        {
            public Arguments()
            {
                this.inputState = "Greedy";
            }

            public string inputState { get; set; }

            public IntegralDataSchema hamiltonian { get; set; }
        }

        public JWEncodeMagic()
        {
            this.Name = $"%jw_encode";
            this.Documentation = new Documentation() { Summary = "Encodes a FermionHamiltonian using he JordanWigner encoding." };
            this.Kind = SymbolKind.Magic;
            this.Execute = this.Run;
        }

        public ExecutionResult Run(string input, IChannel channel)
        {
            var args = Newtonsoft.Json.JsonConvert.DeserializeObject<Arguments>(input);

            var hamiltonian = LoadData.LoadIntegralData(args.hamiltonian).First();

            return JordanWignerEncoding.Create(hamiltonian).QSharpData(args.inputState).ToExecutionResult();
        }
    }
}