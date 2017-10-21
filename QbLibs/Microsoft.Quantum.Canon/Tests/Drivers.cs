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
        public void testDriver()
        {
            var sim = GetSimulator();
            EmptyTest.Run(sim).Wait();
        }
        [TestMethod]
        public void testDriver2()
        {
            var sim = GetSimulator();
            test234.Run(sim).Wait();
        }



        /*
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

    [TestClass]
    public class TestDrivers
    {
        public static IOperationFactory GetSimulator()
        {
            return new QuantumSimulator();
        }

        [TestMethod]
        public void testDriver3()
        {
            var sim = GetSimulator();
            var msg = test234.Run(sim).Result;
            msg = Result.Zero;
            msg = Result.One;
        }



        //var toffoli = new ToffoliSimulator();
        //var result = HelloWorld.Run(toffoli, 12).Result;
        //System.Console.Out.WriteLine(result);

        //var qsim = new QuantumSimulator();
        //var msg = Result.One;
        //test123.Run(qsim).Wait();
        //var r = TeleportTest.Run(qsim, msg).Result;
        //System.Console.WriteLine($"Sent: {msg}; Received: {r}");

        //msg = Result.Zero;
        //r = TeleportTest.Run(qsim, msg).Result;
        //System.Console.WriteLine($"Sent: {msg}; Received: {r}");

        //System.Console.ReadKey();
    }
}
