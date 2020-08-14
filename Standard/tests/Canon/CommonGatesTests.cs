// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.


using Microsoft.Quantum.Simulation.Core;
using Microsoft.Quantum.Simulation.Simulators;
using System.Collections.Generic;
using Xunit;

namespace Microsoft.Quantum.Tests
{
    internal class MockQubit : Qubit
    {
        public MockQubit(int id) : base(id)
        {
        }
    }

    public class CommonGatesRuntimeMetadataTests
    {
        [Fact]
        public void CX()
        {
            var control = new MockQubit(1);
            var target = new MockQubit(0);
            var op = new Microsoft.Quantum.Canon.CX(new QuantumSimulator());
            var args = op.__dataIn((control, target));
            var expected = new RuntimeMetadata()
            {
                Label = "X",
                FormattedNonQubitArgs = "",
                IsAdjoint = false,
                IsControlled = true,
                IsMeasurement = false,
                IsComposite = false,
                Children = null,
                Controls = new List<Qubit>() { control },
                Targets = new List<Qubit>() { target },
            };

            Assert.Equal(op.GetRuntimeMetadata(args), expected);
        }

        [Fact]
        public void CY()
        {
            var control = new MockQubit(1);
            var target = new MockQubit(0);
            var op = new Microsoft.Quantum.Canon.CY(new QuantumSimulator());
            var args = op.__dataIn((control, target));
            var expected = new RuntimeMetadata()
            {
                Label = "Y",
                FormattedNonQubitArgs = "",
                IsAdjoint = false,
                IsControlled = true,
                IsMeasurement = false,
                IsComposite = false,
                Children = null,
                Controls = new List<Qubit>() { control },
                Targets = new List<Qubit>() { target },
            };

            Assert.Equal(op.GetRuntimeMetadata(args), expected);
        }

        [Fact]
        public void CZ()
        {
            var control = new MockQubit(1);
            var target = new MockQubit(0);
            var op = new Microsoft.Quantum.Canon.CZ(new QuantumSimulator());
            var args = op.__dataIn((control, target));
            var expected = new RuntimeMetadata()
            {
                Label = "Z",
                FormattedNonQubitArgs = "",
                IsAdjoint = false,
                IsControlled = true,
                IsMeasurement = false,
                IsComposite = false,
                Children = null,
                Controls = new List<Qubit>() { control },
                Targets = new List<Qubit>() { target },
            };

            Assert.Equal(op.GetRuntimeMetadata(args), expected);
        }
    }
}
