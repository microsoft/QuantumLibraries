// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

using System;

using Microsoft.Quantum.Simulation.Common;
using Microsoft.Quantum.Simulation.Simulators;
using Xunit;
using Assert = Xunit.Assert;

namespace Microsoft.Quantum.Tests
{
    public class EstimateFrequency
    {
        [Fact]
        public void EstimateFrequencyWithAndWithoutEmulation()
        {
            void TestOne(SimulatorBase sim, int expected)
            {
                var count = 0;
                sim.OnLog += (_) => count++;
                sim.DisableLogToConsole();

                EstimateFrequencyEmulationTest.Run(sim).Wait();
                Assert.Equal(expected, count);

                if (sim is IDisposable dis) dis.Dispose();
            }

            TestOne(new QuantumSimulator(), 1);
            TestOne(new ToffoliSimulator(), 2000);
        }

        [Fact]
        public void TestEstimateFrequencyBinomial()
        {
            using var sim = new QuantumSimulator(randomNumberGeneratorSeed: 655321);
            TestEstimateFrequencyBinomialInner.Run(sim).Wait();
        }

        [Fact]
        public void TestRobustPhaseEstimation()
        {
            using var sim = new QuantumSimulator(randomNumberGeneratorSeed: 655321);
            TestRobustPhaseEstimationInner.Run(sim).Wait();
        }
    }
}
