using System;
using Microsoft.Quantum.Simulation.Simulators;
using Microsoft.Quantum.Simulation.Core;
using System.Threading.Tasks;
using System.Collections.Generic;
using System.Reflection;
using System.Linq;
using Xunit;
using Xunit.Sdk;
using Xunit.Abstractions;
using System.Runtime.Serialization;
using System.ComponentModel;
using System.Threading;

namespace Microsoft.Quantum.Canon.Tests
{
    public class Drivers
    {

        private readonly ITestOutputHelper TestOutputHelper;
        public Drivers(ITestOutputHelper testOutputHelper)
        {
            TestOutputHelper = testOutputHelper;
        }



        public static IOperationFactory GetSimulator()
        {
            return new QuantumSimulator();
        }
        
        [Theory(DisplayName ="OperationTest")]
        [OperationData(suffix: "Test")]
        async public Task OperationDriver(String name, TypeInfo operation)
        {
            var sim = GetSimulator();
            await (typeof(QuantumSimulator)
                       .GetMethod("Run")
                       .MakeGenericMethod(operation, typeof(QVoid), typeof(QVoid))
                       .Invoke(sim, new[] { QVoid.Instance })
                   as Task);
        }

        [Theory(DisplayName = "EXPECTED FAIL: OperationTestShouldFail")]
        [OperationData(suffix: "TestExFail")]
        async public Task OperationExFailDriver(String name, TypeInfo operation)
        {
            var sim = GetSimulator();
            var threw = false;
            try
            {
                await (typeof(QuantumSimulator)
                           .GetMethod("Run")
                           .MakeGenericMethod(operation, typeof(QVoid), typeof(QVoid))
                           .Invoke(sim, new[] { QVoid.Instance })
                       as Task);
            }
            catch (ExecutionFailException ex)
            {
                threw = true;
                TestOutputHelper.WriteLine($"EXPECTED FAILURE: {ex.Message}");
            }

            if (!threw)
            {
                throw new Exception($"TEST FAILURE RESOLVED. Rename test to remove ExFail.");
            }
        }

        [Theory(DisplayName = "OperationTestShouldFail")]
        [OperationData(suffix: "TestShouldFail")]
        async public Task OperationShouldFailDriver(String name, TypeInfo operation)
        {
            var sim = GetSimulator();
            var threw = false;
            try
            {
                await (typeof(QuantumSimulator)
                           .GetMethod("Run")
                           .MakeGenericMethod(operation, typeof(QVoid), typeof(QVoid))
                           .Invoke(sim, new[] { QVoid.Instance })
                       as Task);
            } catch (ExecutionFailException ex)
            {
                threw = true;
            }

            if (!threw)
            {
                throw new Exception("Operation did not fail, but should have.");
            }
        }

    }
}
