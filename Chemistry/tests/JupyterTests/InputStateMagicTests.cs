// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using System.Collections.Generic;
using System.IO;
using System.Linq;
using Microsoft.Jupyter.Core;
using Microsoft.Quantum.Chemistry.Broombridge;
using Microsoft.Quantum.Chemistry.Fermion;
using Microsoft.Quantum.Chemistry.Magic;
using Newtonsoft.Json;
using Microsoft.Quantum.Chemistry.LadderOperators;
using Xunit;
using Microsoft.Quantum.Chemistry.Generic;

namespace Microsoft.Quantum.Chemistry.Tests
{
    public class InputStateMagicTests
    {
        public (InputStateMagic, MockChannel) Init() =>
            (new InputStateMagic(), new MockChannel());

        [Fact]
        public void LoadNoInput()
        {
            var (magic, channel) = Init();

            Assert.Equal("%inputstate-load", magic.Name);
            var result = magic.Run("", channel);

            Assert.Equal(ExecuteStatus.Error, result.Status);
        }

        [Fact]
        public void LoadInvalidFile()
        {
            var (magic, channel) = Init();
            var args = JsonConvert.SerializeObject(new InputStateMagic.Arguments
            {
                fileName = "foo_bar.yaml"
            });

            Assert.Throws<FileNotFoundException>(() => magic.Run(args, channel));
        }

        [Fact]
        public void LoadFromBroombridgeFile()
        {
            var (magic, channel) = Init();
            var args = JsonConvert.SerializeObject(new InputStateMagic.Arguments
            {
                fileName = "broombridge_v0.2.yaml",
                wavefunctionLabel = "UCCSD |G>"
            });

            var result = magic.Run(args, channel);
            var wavefunction = (InputState)result.Output;
            Assert.Equal(ExecuteStatus.Ok, result.Status);
            Assert.Equal("UCCSD |G>", wavefunction.Label);
            Assert.Equal(StateType.UnitaryCoupledCluster, wavefunction.TypeOfState);
        }


        [Fact]
        public void LoadFromProblemDescription()
        {
            var (magic, channel) = Init();
            var broombridgeMagic = new BroombridgeMagic();
            var broombridge = (CurrentVersion.Data)broombridgeMagic.Run("broombridge_v0.2.yaml", channel).Output;

            var args = JsonConvert.SerializeObject(new FermionHamiltonianLoadMagic.Arguments
            {
                problemDescription = broombridge.ProblemDescriptions.First()
            });

            var result = magic.Run(args, channel);
            var wavefunction = (InputState)result.Output;
            Assert.Equal(ExecuteStatus.Ok, result.Status);
            Assert.Equal("Greedy", wavefunction.Label);
            Assert.Equal(StateType.SparseMultiConfigurational, wavefunction.TypeOfState);
        }
    }
}