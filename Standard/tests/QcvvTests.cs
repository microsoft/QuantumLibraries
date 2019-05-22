using System;
using System.Linq;
using System.Runtime.InteropServices;

using Microsoft.Quantum.Intrinsic;
using Microsoft.Quantum.Simulation;
using Microsoft.Quantum.Simulation.Common;
using Microsoft.Quantum.Simulation.Core;
using Microsoft.Quantum.Simulation.Simulators;
using Microsoft.Quantum.Standard.Emulation;
using Xunit;
using Assert = Xunit.Assert;

namespace Microsoft.Quantum.Tests
{
    public class EstimateFrequency
    {
        [Fact]
        public void TestEmulation()
        {
            void TestOne(SimulatorBase sim, int expected)
            {
                var count = 0;
                sim.OnLog += (_) => count++;

                EstimateFrequencyEmulationTest.Run(sim).Wait();
                Assert.Equal(expected, count);

                if (sim is IDisposable dis) dis.Dispose();
            }

            TestOne(new QuantumSimulator(), 1);
            TestOne(new ToffoliSimulator(), 2000);
        }
    }
}
