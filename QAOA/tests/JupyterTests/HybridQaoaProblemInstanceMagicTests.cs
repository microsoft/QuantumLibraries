// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using System.Threading.Tasks;
using Microsoft.Jupyter.Core;
using Newtonsoft.Json;
using Xunit;
using System;
using Microsoft.Quantum.QAOA.JupyterTests;
using Microsoft.Quantum.QAOA.QaoaHybrid;
using QAOA.Jupyter;

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

            double[] dh = new Double[] { 0, 0 };
            double[] dJ = new Double[]{ 0, 1,
                               0, 0};

            var args = JsonConvert.SerializeObject(new HybridQaoaProblemInstanceMagic.Arguments
            {
                OneLocalHamiltonianCoefficients = dh,
                TwoLocalHamiltonianCoefficients = dJ
            });

            var result = await magic.Run(args, channel);
            var problemInstance = result.Output as ProblemInstance;
            Assert.Equal(ExecuteStatus.Ok, result.Status);

            Assert.Equal(problemInstance.problemSizeInBits, dh.Length);
            Assert.Equal(problemInstance.oneLocalHamiltonianCoefficients, dh);
            Assert.Equal(problemInstance.twoLocalHamiltonianCoefficients,dJ);
        }

    }
}