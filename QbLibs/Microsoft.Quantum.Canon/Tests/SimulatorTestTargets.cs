using Microsoft.Quantum.Testing;
using Microsoft.Quantum.Canon.Tests;
using Microsoft.Quantum.Simulation.Core;
using Xunit;

namespace Microsoft.Quantum.Canon
{
    public class SimulatorTestTargets
    {
        // note that one can provide custom namespace where to search for tests
        [QuantumTestTarget(TestNamespace = "Microsoft.Quantum.Canon")]
        public void QuantumSimulatorTarget( QuantumTestCaseOperationData opData )
        {
            var sim = TestBase.GetSimulator();
            opData.OperationRunner(sim);
        }

        // when Skip attribute is specified the test will be skipped and marked with yellow icon
        [QuantumTestTarget(Suffix = "TestExFail", Skip ="This is the test that is known to fail")]
        public void QuantumSimulatorTargetExFail(QuantumTestCaseOperationData opData)
        {
            var sim = TestBase.GetSimulator();
            opData.OperationRunner(sim);
        }

        // one can find tests with custom suffix
        [QuantumTestTarget(Suffix = "TestShouldFail")]
        public void QuantumSimulatorTargetShouldFail(QuantumTestCaseOperationData opData)
        {
            var sim = TestBase.GetSimulator();
            Assert.ThrowsAny<ExecutionFailException>(() => opData.OperationRunner(sim));
        }
    }
}
