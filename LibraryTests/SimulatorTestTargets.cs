// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using Microsoft.Quantum.Testing;
using Microsoft.Quantum.Simulation.Core;
using Xunit;
using Microsoft.Quantum.Simulation.Simulators;

namespace Microsoft.Quantum.Tests
{
    public class SimulatorTestTargets
    {
        // note that one can provide custom namespace where to search for tests
        [QuantumTestTarget(TestNamespace = "Microsoft.Quantum.Tests")]
        public void QuantumSimulatorTarget(QuantumTestCaseOperationData opData)
        {
            using (var sim = new QuantumSimulator())
            {
                opData.OperationRunner(sim);
            }
        }

        [QuantumTestTarget(TestNamespace = "Microsoft.Quantum.Canon")]
        public void QuantumSimulatorCanonTarget(QuantumTestCaseOperationData opData)
        {
            using (var sim = new QuantumSimulator())
            {
                opData.OperationRunner(sim);
            }
        }

        [QuantumTestTarget(TestNamespace = "Microsoft.Quantum.Canon", AssemblyName ="Microsoft.Quantum.Canon")]
        public void QuantumSimulatorOldCanonTarget(QuantumTestCaseOperationData opData)
        {
            using (var sim = new QuantumSimulator())
            {
#if DEBUG
                opData.OperationRunner(sim);
#else
                throw new System.Exception("Tests should be removed from Microsoft.Quantum.Canon.dll.");
#endif
            }
        }

        // when Skip attribute is specified the test will be skipped and marked with yellow icon
        [QuantumTestTarget(Suffix = "TestExFail", Skip = "This test is expected to fail.")]
        public void QuantumSimulatorTargetExFail(QuantumTestCaseOperationData opData)
        {
            using (var sim = new QuantumSimulator())
            {
                opData.OperationRunner(sim);
            }
        }

        // one can find tests with custom suffix
        [QuantumTestTarget(Suffix = "TestShouldFail")]
        public void QuantumSimulatorTargetShouldFail(QuantumTestCaseOperationData opData)
        {
            using (var sim = new QuantumSimulator())
            {
                Assert.ThrowsAny<ExecutionFailException>(() => opData.OperationRunner(sim));
            }
        }

        // FIXME: there has to be a cleaner way of doing this.
        [QuantumTestTarget(TestNamespace = "Microsoft.Quantum.Samples.DatabaseSearch", AssemblyName ="DatabaseSearchSample")]
        public void DatabaseSearchTarget(QuantumTestCaseOperationData opData)
        {
            using (var sim = new QuantumSimulator())
            {
                opData.OperationRunner(sim);
            }
        }

    }
}
