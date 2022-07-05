// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.


using Microsoft.Quantum.Simulation.Simulators;
using Xunit;

namespace Microsoft.Quantum.Tests
{
    public class MultiplexerTests
    {
        [Fact]
        public void TestMultiplexerResources()
        {
            foreach (var (numControls, expectedTCount, expectedWidth) in new [] {
                    (1, 0, 6),
                    (2, 8, 8),
                    (3, 24, 10),
                    (4, 56, 12),
                    (5, 120, 14),
                    (6, 248, 16)
                }) {
                var estimator = new ResourcesEstimator();
                EstimateMultiplexOperationsCosts.Run(estimator, numControls, 1 << numControls).Wait();
                var tcount = (double)estimator.Data.Rows.Find("T")[1];
                var width = (double)estimator.Data.Rows.Find("Width")[1];
                Assert.Equal(expectedTCount, tcount);
                Assert.Equal(expectedWidth, width);
            }
        }
    }
}
