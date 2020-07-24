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
    public class ForEachRuntimeMetadataTests
    {
        [Fact]
        public void ForEach()
        {
            var op = new Microsoft.Quantum.Arrays.ForEach<Qubit, Result>(new QuantumSimulator());
            var baseOp = new Microsoft.Quantum.Intrinsic.M(new QuantumSimulator());
            IQArray<Qubit> targets = new QArray<Qubit>(new Qubit[] { });
            var args = op.__dataIn((baseOp, targets));
            var expected = new RuntimeMetadata()
            {
                Label = "ForEach",
                FormattedNonQubitArgs = "(M)",
                IsComposite = true,
                Targets = new List<Qubit>() { },
            };
            Assert.Equal(expected, op.GetRuntimeMetadata(args));
        }
    }
}
