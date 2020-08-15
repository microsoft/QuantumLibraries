// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using System.Threading.Tasks;
using Microsoft.Jupyter.Core;
using Newtonsoft.Json;
using System;
using Microsoft.Quantum.QAOA.QaoaHybrid;
using Xunit;
using Microsoft.Quantum.QAOA.Jupyter;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Assert = Xunit.Assert;

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
            var oneLocalHamiltonianCoefficients = new double[] { 0, 0 };
            var twoLocalHamiltonianCoefficients = new double[] { 0, 1, 0, 0};
            var initialBeta = new double[] {0, 0};
            var initialGamma = new double[] { 0, 0 };

            var simpleMaxCut = new ProblemInstance(oneLocalHamiltonianCoefficients, twoLocalHamiltonianCoefficients);

            var args = JsonConvert.SerializeObject(new HybridQaoaRunMagic.Arguments
            {
                NumberOfIterations = numberOfIterations,
                p = p,
                ProblemInstance = simpleMaxCut,
                InitialBeta = initialBeta,
                InitialGamma = initialGamma
            });

            var result = await magic.Run(args, channel);
            var optimalSolution = result.Output as Solution;
            Assert.Equal(ExecuteStatus.Ok, result.Status);
            var optimizationResult1 = new[] {false, true};
            var optimizationResult2 = new[] {true, false};

            if (optimalSolution.SolutionVector[0] == false)
            {
                CollectionAssert.AreEqual(optimalSolution.SolutionVector, optimizationResult1, "Hybrid QAOA produced incorrect result when running magic.");
            }
            else
            {
                CollectionAssert.AreEqual(optimalSolution.SolutionVector, optimizationResult2, "Hybrid QAOA produced incorrect result when running magic.");
            }
        }
    }

    public class HybridQaoaWithRandomParametersMagicTests
    {
        public (HybridQaoaWithRandomParametersRunMagic, MockChannel) Init() =>
            (new HybridQaoaWithRandomParametersRunMagic(), new MockChannel());

        [Fact]
        public async Task HybridQaoaRun()
        {
            var (magic, channel) = Init();
            Assert.Equal("%qaoa.hybridqaoa.random.params.run", magic.Name);

            var numberOfIterations = 50;
            var p = 2;
            var oneLocalHamiltonianCoefficients = new double[] { 0, 0 };
            var twoLocalHamiltonianCoefficients = new double[] { 0, 1, 0, 0 };
            var numberOfRandomStartingPoints = 2;

            var simpleMaxCut = new ProblemInstance(oneLocalHamiltonianCoefficients, twoLocalHamiltonianCoefficients);

            var args = JsonConvert.SerializeObject(new HybridQaoaWithRandomParametersRunMagic.Arguments
            {
                NumberOfIterations = numberOfIterations,
                p = p,
                ProblemInstance = simpleMaxCut,
                NumberOfRandomStartingPoints = numberOfRandomStartingPoints
        });

            var result = await magic.Run(args, channel);
            var optimalSolution = result.Output as Solution;
            Assert.Equal(ExecuteStatus.Ok, result.Status);
            var optimizationResult1 = new[] { false, true };
            var optimizationResult2 = new[] { true, false };

            if (optimalSolution.SolutionVector[0] == false)
            {
                CollectionAssert.AreEqual(optimalSolution.SolutionVector, optimizationResult1, "Hybrid QAOA produced incorrect result when running magic.");
            }
            else
            {
                CollectionAssert.AreEqual(optimalSolution.SolutionVector, optimizationResult2, "Hybrid QAOA produced incorrect result when running magic.");
            }
        }
    }

}