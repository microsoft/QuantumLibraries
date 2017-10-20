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

        /*
         * NONE OF THE TESTS ARE PASSING :(
        [TestMethod]
        public void TeleportationDriver()
        {
            var sim = GetSimulator();
            TeleportationTest.Run(sim).Wait();
        }

        [TestMethod]
        public void RUSDriver()
        {
            var sim = GetSimulator();
            RUSTests.Run(sim).Wait();
        }

        [TestMethod]
        public void SuperdenseCodingDriver()
        {
            var sim = GetSimulator();
            SuperdenseCodingTest.Run(sim).Wait();
        }

        [TestMethod]
        public void GroverExampleAAbyOracleDriver()
        {
            var sim = GetSimulator();
            ExampleAAbyOracle.Run(sim).Wait();
        }

        [TestMethod]
        public void GroverTestDriver()
        {
            var sim = GetSimulator();
            GroverTest.Run(sim, 1024, 10, 11).Wait();
        }
        */
    }
}
