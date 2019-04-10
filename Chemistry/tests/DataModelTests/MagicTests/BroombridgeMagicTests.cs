// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using System.Collections.Generic;
using System.IO;
using System.Linq;
using Microsoft.Jupyter.Core;
using Microsoft.Quantum.Chemistry.Broombridge;
using Microsoft.Quantum.Chemistry.Magic;
using Newtonsoft.Json;

using Xunit;

namespace Microsoft.Quantum.Chemistry.Tests
{
    public class BroombridgeMagicTests
    {
        [Fact]
        public void LoadNoFile()
        {
            var magic = new BroombridgeMagic();
            var channel = new MockChannel();

            Assert.Equal("%broombridge", magic.Name);
            var result = magic.Run("", channel);

            Assert.Equal(ExecuteStatus.Error, result.Status);
        }


        [Fact]
        public void LoadInvalidFile()
        {
            var filename = "Broombridge/foo_bar.yaml";
            var magic = new BroombridgeMagic();
            var channel = new MockChannel();

            Assert.Throws<FileNotFoundException>(() => magic.Run(filename, channel));
        }


        [Fact]
        public void LoadBroombridgeFile()
        {
            var filename = "Broombridge/broombridge_v0.2.yaml";
            var magic = new BroombridgeMagic();
            var channel = new MockChannel();

            var result = magic.Run(filename, channel);
            var broombridge = (CurrentVersion.Data)result.Output;
            Assert.Equal(ExecuteStatus.Ok, result.Status);
            Assert.Equal("0.2", broombridge.Format.Version);
            Assert.Equal(3, broombridge.Bibliography.Count);
            Assert.Equal(1, broombridge.ProblemDescriptions.Count);
        }
    }
}