// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Qaoa.JupyterTests
{
    using System.Threading.Tasks;
    using Microsoft.Jupyter.Core;
    using Newtonsoft.Json;
    using Xunit;
    using Microsoft.Quantum.Qaoa.QaoaHybrid;
    using Microsoft.Quantum.Qaoa.Jupyter;

    public class HybridQaoaParametersMagicTests
    {
        public (HybridQaoaParametersMagic, MockChannel) Init() =>
            (new HybridQaoaParametersMagic(), new MockChannel());

        [Fact]
        public async Task HybridQaoaProblemInstance()
        {
            var (magic, channel) = Init();
            Assert.Equal("%qaoa.hybridqaoa.create.parameters", magic.Name);

            var beta = new double[] { 1, 2, 3 };
            var gamma = new double[] { 5, 6, 7 };


            var args = JsonConvert.SerializeObject(new HybridQaoaParametersMagic.Arguments
            {
                Beta = beta,
                Gamma = gamma
            });

            var result = await magic.Run(args, channel);
            var qaoaParameters = result.Output as QaoaParameters;
            Assert.Equal(ExecuteStatus.Ok, result.Status);

            Assert.Equal(qaoaParameters.Betas, beta);
            Assert.Equal(qaoaParameters.Gammas, gamma);
        }

    }
}