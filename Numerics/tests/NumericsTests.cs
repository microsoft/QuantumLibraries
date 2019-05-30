// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using Microsoft.Quantum.Simulation.XUnit;
using Microsoft.Quantum.Simulation.Simulators;
using Xunit.Abstractions;
using System.Diagnostics;

namespace Microsoft.Quantum.Numerics.Tests
{
    public class NumericsTests
    {
        private readonly ITestOutputHelper output;

        public NumericsTests(ITestOutputHelper output)
        {
            this.output = output;
        }

        /// <summary>
        /// This driver will run all Q# tests (operations named "...Test")
        /// that are located inside Microsoft.Quantum.Numerics.Tests using the quantum
        /// simulator.
        /// </summary>
        [OperationDriver(TestNamespace = "Microsoft.Quantum.Numerics.Tests",
                         TestCasePrefix = "QSim:")]
        public void QSimTests(TestOperation op)
        {
            var sim = new QuantumSimulator();
            // OnLog defines action(s) performed when Q# test calls function Message
            sim.OnLog += (msg) => { output.WriteLine(msg); };
            sim.OnLog += (msg) => { Debug.WriteLine(msg); };
            op.TestOperationRunner(sim);
        }

        /// <summary>
        /// This driver will run all Q# tests (operations named "...Test")
        /// that are located inside Microsoft.Quantum.Numerics.ToffoliTests using the
        /// Toffoli simulator.
        /// </summary>
        [OperationDriver(TestNamespace = "Microsoft.Quantum.Numerics.ToffoliTests",
                         TestCasePrefix = "ToffSim:")]
        public void ToffoliSimTests(TestOperation op)
        {
            var sim = new ToffoliSimulator();
            // OnLog defines action(s) performed when Q# test calls function Message
            sim.OnLog += (msg) => { output.WriteLine(msg); };
            sim.OnLog += (msg) => { Debug.WriteLine(msg); };
            op.TestOperationRunner(sim);
        }
    }
}
