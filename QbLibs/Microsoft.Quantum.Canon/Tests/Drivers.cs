using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Microsoft.Quantum.Simulation.Simulators;
using Microsoft.Quantum.Simulation.Core;

namespace Microsoft.Quantum.Canon.Tests
{
    [TestClass]
    public class Drivers
    {
        public static IOperationFactory GetSimulator()
        {
            return new QuantumSimulator();
        }

        [TestMethod]
        public void TeleportationDriver()
        {
            var sim = GetSimulator();
            TeleportationTest.Run(sim);
        }

        [TestMethod]
        public void RUSDriver()
        {
            var sim = GetSimulator();
            TeleportationTest.Run(sim);
        }

        [TestMethod]
        public void SuperdenseCodingDriver()
        {
            var sim = GetSimulator();
            SuperdenseCodingTest.Run(sim);
        }
    }
}
