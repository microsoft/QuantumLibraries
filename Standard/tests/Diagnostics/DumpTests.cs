// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.InteropServices;
using Microsoft.Quantum.Diagnostics.Emulation;
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
    public class DumpOperation
    {
        [Fact]
        public void DumpSReturnsCorrectMatrix()
        {
            using var sim = new QuantumSimulator();
            var diagnostics = new List<object>();
            sim.OnDisplayableDiagnostic += diagnostic =>
            {
                diagnostics.Add(diagnostic);
            };
            sim.DisableLogToConsole();

            DumpS.Run(sim).Wait();

            Assert.Equal(diagnostics.Count, 1);
            var diagnostic = diagnostics.Single();

            Assert.IsType<DisplayableUnitaryOperator>(diagnostic);
            var unitary = diagnostic as DisplayableUnitaryOperator;

            Assert.NotNull(unitary);
            Assert.Equal(1, unitary.Qubits?.Count);
            Assert.NotNull(unitary.Data);
            Assert.Equal(3, unitary.Data.ndim);
            Assert.Equal(2, unitary.Data.shape[0]);
            Assert.Equal(2, unitary.Data.shape[1]);
            Assert.Equal(2, unitary.Data.shape[2]);
            // Check that the matrix is [1 0; 0 𝑖].
            Assert.Equal(1.0, (double)unitary.Data[0, 0, 0], precision: 6);
            Assert.Equal(0.0, (double)unitary.Data[0, 1, 0], precision: 6);
            Assert.Equal(0.0, (double)unitary.Data[1, 0, 0], precision: 6);
            Assert.Equal(0.0, (double)unitary.Data[1, 1, 0], precision: 6);
            Assert.Equal(0.0, (double)unitary.Data[0, 0, 1], precision: 6);
            Assert.Equal(0.0, (double)unitary.Data[0, 1, 1], precision: 6);
            Assert.Equal(0.0, (double)unitary.Data[1, 0, 1], precision: 6);
            Assert.Equal(1.0, (double)unitary.Data[1, 1, 1], precision: 6);
        }

        
        [Fact]
        public void DumpCnotReturnsCorrectMatrix()
        {
            using var sim = new QuantumSimulator();
            var diagnostics = new List<object>();
            sim.OnDisplayableDiagnostic += diagnostic =>
            {
                diagnostics.Add(diagnostic);
            };
            sim.DisableLogToConsole();

            DumpCnot.Run(sim).Wait();

            Assert.Equal(diagnostics.Count, 1);
            var diagnostic = diagnostics.Single();

            Assert.IsType<DisplayableUnitaryOperator>(diagnostic);
            var unitary = diagnostic as DisplayableUnitaryOperator;

            Assert.NotNull(unitary);
            Assert.Equal(2, unitary.Qubits?.Count);
            Assert.NotNull(unitary.Data);
            Assert.Equal(3, unitary.Data.ndim);
            Assert.Equal(4, unitary.Data.shape[0]);
            Assert.Equal(4, unitary.Data.shape[1]);
            Assert.Equal(2, unitary.Data.shape[2]);
            // Check that the matrix is [1 0 0 0; 0 0 0 1; 0 0 1 0; 0 1 0 0].
            // real
            Assert.Equal(1.0, (double)unitary.Data[0, 0, 0], precision: 6);
            Assert.Equal(0.0, (double)unitary.Data[0, 1, 0], precision: 6);
            Assert.Equal(0.0, (double)unitary.Data[0, 2, 0], precision: 6);
            Assert.Equal(0.0, (double)unitary.Data[0, 3, 0], precision: 6);
            Assert.Equal(0.0, (double)unitary.Data[1, 0, 0], precision: 6);
            Assert.Equal(0.0, (double)unitary.Data[1, 1, 0], precision: 6);
            Assert.Equal(0.0, (double)unitary.Data[1, 2, 0], precision: 6);
            Assert.Equal(1.0, (double)unitary.Data[1, 3, 0], precision: 6);
            Assert.Equal(0.0, (double)unitary.Data[2, 0, 0], precision: 6);
            Assert.Equal(0.0, (double)unitary.Data[2, 1, 0], precision: 6);
            Assert.Equal(1.0, (double)unitary.Data[2, 2, 0], precision: 6);
            Assert.Equal(0.0, (double)unitary.Data[2, 3, 0], precision: 6);
            Assert.Equal(0.0, (double)unitary.Data[3, 0, 0], precision: 6);
            Assert.Equal(1.0, (double)unitary.Data[3, 1, 0], precision: 6);
            Assert.Equal(0.0, (double)unitary.Data[3, 2, 0], precision: 6);
            Assert.Equal(0.0, (double)unitary.Data[3, 3, 0], precision: 6);

            // imaginary
            Assert.Equal(0.0, (double)unitary.Data[0, 0, 1], precision: 6);
            Assert.Equal(0.0, (double)unitary.Data[0, 1, 1], precision: 6);
            Assert.Equal(0.0, (double)unitary.Data[0, 2, 1], precision: 6);
            Assert.Equal(0.0, (double)unitary.Data[0, 3, 1], precision: 6);
            Assert.Equal(0.0, (double)unitary.Data[1, 0, 1], precision: 6);
            Assert.Equal(0.0, (double)unitary.Data[1, 1, 1], precision: 6);
            Assert.Equal(0.0, (double)unitary.Data[1, 2, 1], precision: 6);
            Assert.Equal(0.0, (double)unitary.Data[1, 3, 1], precision: 6);
            Assert.Equal(0.0, (double)unitary.Data[2, 0, 1], precision: 6);
            Assert.Equal(0.0, (double)unitary.Data[2, 1, 1], precision: 6);
            Assert.Equal(0.0, (double)unitary.Data[2, 2, 1], precision: 6);
            Assert.Equal(0.0, (double)unitary.Data[2, 3, 1], precision: 6);
            Assert.Equal(0.0, (double)unitary.Data[3, 0, 1], precision: 6);
            Assert.Equal(0.0, (double)unitary.Data[3, 1, 1], precision: 6);
            Assert.Equal(0.0, (double)unitary.Data[3, 2, 1], precision: 6);
            Assert.Equal(0.0, (double)unitary.Data[3, 3, 1], precision: 6);
        }
    }
}
