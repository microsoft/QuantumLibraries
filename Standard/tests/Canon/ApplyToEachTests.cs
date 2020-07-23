// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.


using Microsoft.Quantum.Simulation.Core;
using Microsoft.Quantum.Simulation.Simulators;
using Microsoft.Quantum.Simulation.XUnit;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Security.Cryptography;
using System.Text;
using Xunit;

namespace Microsoft.Quantum.Tests
{
    public class ApplyToEachRuntimeMetadataTests
    {
        [Fact]
        public void ApplyToEach()
        {
            var op = new Microsoft.Quantum.Canon.ApplyToEach<Qubit>(new QuantumSimulator());
            var baseOp = new Microsoft.Quantum.Intrinsic.I(new QuantumSimulator());
            IQArray<Qubit> targets = new QArray<Qubit>(new Qubit[] { });
            var args = op.__dataIn((baseOp, targets));
            var expected = new RuntimeMetadata() { IsComposite = true };
            Assert.Equal(expected, op.GetRuntimeMetadata(args));
        }

        [Fact]
        public void ApplyToEachC()
        {
            var op = new Microsoft.Quantum.Canon.ApplyToEachC<Qubit>(new QuantumSimulator());
            var baseOp = new Microsoft.Quantum.Intrinsic.I(new QuantumSimulator());
            IQArray<Qubit> targets = new QArray<Qubit>(new Qubit[] { });
            var args = op.__dataIn((baseOp, targets));
            var expected = new RuntimeMetadata() { IsComposite = true };
            Assert.Equal(expected, op.GetRuntimeMetadata(args));
        }

        [Fact]
        public void ApplyToEachA()
        {
            var op = new Microsoft.Quantum.Canon.ApplyToEachA<Qubit>(new QuantumSimulator());
            var baseOp = new Microsoft.Quantum.Intrinsic.I(new QuantumSimulator());
            IQArray<Qubit> targets = new QArray<Qubit>(new Qubit[] { });
            var args = op.__dataIn((baseOp, targets));
            var expected = new RuntimeMetadata() { IsComposite = true };
            Assert.Equal(expected, op.GetRuntimeMetadata(args));
        }

        [Fact]
        public void ApplyToEachCA()
        {
            var op = new Microsoft.Quantum.Canon.ApplyToEachCA<Qubit>(new QuantumSimulator());
            var baseOp = new Microsoft.Quantum.Intrinsic.I(new QuantumSimulator());
            IQArray<Qubit> targets = new QArray<Qubit>(new Qubit[] { });
            var args = op.__dataIn((baseOp, targets));
            var expected = new RuntimeMetadata() { IsComposite = true };
            Assert.Equal(expected, op.GetRuntimeMetadata(args));
        }
    }
}
