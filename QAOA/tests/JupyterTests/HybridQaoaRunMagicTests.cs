// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using System.Threading.Tasks;
using Microsoft.Jupyter.Core;
using Newtonsoft.Json;
using System;
using Microsoft.Quantum.QAOA.QaoaHybrid;
using Xunit;
using Microsoft.Quantum.QAOA.Jupyter;

namespace Microsoft.Quantum.QAOA.JupyterTests
{
    public class HybridQaoaRunMagicTests
    {
        public (HybridQaoaRunMagic, MockChannel) Init() =>
            (new HybridQaoaRunMagic(), new MockChannel());

        [Fact]
        public async Task HybridQaoaRun()
        {
            var (magic, channel) = Init();
            Assert.Equal("%qaoa.hybridqaoa.run", magic.Name);

            var numberOfIterations = 50;
            var p = 2;
            double[] oneLocalHamiltonianCoefficients = new Double[] { 0, 0 };
            double[] twoLocalHamiltonianCoefficients = new Double[]{ 0, 1, 0, 0};

            ProblemInstance simpleMaxCut = new ProblemInstance(oneLocalHamiltonianCoefficients, twoLocalHamiltonianCoefficients);

            var args = JsonConvert.SerializeObject(new HybridQaoaRunMagic.Arguments
            {
                NumberOfIterations = numberOfIterations,
                p = p,
                ProblemInstance = simpleMaxCut
            });

            var result = await magic.Run(args, channel);
            var optimalSolution = result.Output as OptimalSolution;
            Assert.Equal(ExecuteStatus.Ok, result.Status);
            string optimizationResult1 = "01";
            string optimizationResult2 = "10";

            Assert.True(optimalSolution.OptimalVector.Equals(optimizationResult1) || optimalSolution.OptimalVector.Equals(optimizationResult2), "Hybrid QAOA produced incorrect result.");
        }
    }
        
    }