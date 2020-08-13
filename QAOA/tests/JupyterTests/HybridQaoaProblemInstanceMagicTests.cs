// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using System.Threading.Tasks;
using Microsoft.Jupyter.Core;
using Newtonsoft.Json;
using Xunit;
using System;
using Microsoft.Quantum.QAOA.JupyterTests;
using Microsoft.Quantum.QAOA.QaoaHybrid;
using Microsoft.Quantum.QAOA.Jupyter;

namespace Microsoft.Quantum.QAOA.HybridQaoaTests
{

    public class HybridQaoaProblemInstanceMagicTests
    {
        public (HybridQaoaProblemInstanceMagic, MockChannel) Init() =>
        (new HybridQaoaProblemInstanceMagic(), new MockChannel());

        [Fact]
        public async Task HybridQaoaProblemInstance()
        {
            var (magic, channel) = Init();
            Assert.Equal("%qaoa.hybridqaoa.create.problem.instance", magic.Name);

            var oneLocalHamiltonianCoefficients = new double[] {0, 0};
            var twoLocalHamiltonianCoefficients = new double[] {0, 1, 0, 0};

            var args = JsonConvert.SerializeObject(new HybridQaoaProblemInstanceMagic.Arguments
            {
                OneLocalHamiltonianCoefficients = oneLocalHamiltonianCoefficients,
                TwoLocalHamiltonianCoefficients = twoLocalHamiltonianCoefficients
            });

            var result = await magic.Run(args, channel);
            var problemInstance = result.Output as ProblemInstance;
            Assert.Equal(ExecuteStatus.Ok, result.Status);

            Assert.Equal(problemInstance.ProblemSizeInBits, oneLocalHamiltonianCoefficients.Length);
            Assert.Equal(problemInstance.OneLocalHamiltonianCoefficients, oneLocalHamiltonianCoefficients);
            Assert.Equal(problemInstance.TwoLocalHamiltonianCoefficients, twoLocalHamiltonianCoefficients);
        }

    }
}