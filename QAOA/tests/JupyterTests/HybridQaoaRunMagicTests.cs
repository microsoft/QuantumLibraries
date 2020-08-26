// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Qaoa.JupyterTests
{
    using System.Threading.Tasks;
    using Microsoft.Jupyter.Core;
    using Newtonsoft.Json;
    using Microsoft.Quantum.Qaoa.QaoaHybrid;
    using Xunit;
    using Microsoft.Quantum.Qaoa.Jupyter;
    using Microsoft.VisualStudio.TestTools.UnitTesting;
    using Assert = Xunit.Assert;

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
            var initialQaoaParameters = new QaoaParameters(initialBeta, initialGamma);

            var simpleMaxCut = new QaoaProblemInstance(oneLocalHamiltonianCoefficients, twoLocalHamiltonianCoefficients);

            var args = JsonConvert.SerializeObject(new HybridQaoaRunMagic.Arguments
            {
                NumberOfIterations = numberOfIterations,
                p = p,
                QaoaProblemInstance = simpleMaxCut,
                InitialQaoaParameters = initialQaoaParameters,
            });

            var result = await magic.Run(args, channel);
            var optimalSolution = result.Output as QaoaSolution;
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

            var simpleMaxCut = new QaoaProblemInstance(oneLocalHamiltonianCoefficients, twoLocalHamiltonianCoefficients);

            var args = JsonConvert.SerializeObject(new HybridQaoaWithRandomParametersRunMagic.Arguments
            {
                NumberOfIterations = numberOfIterations,
                p = p,
                QaoaProblemInstance = simpleMaxCut,
                NumberOfRandomStartingPoints = numberOfRandomStartingPoints
        });

            var result = await magic.Run(args, channel);
            var optimalSolution = result.Output as QaoaSolution;
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