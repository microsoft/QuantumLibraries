﻿// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Threading.Tasks;
using Microsoft.Jupyter.Core;
using Microsoft.Quantum.Chemistry.Broombridge;
using Microsoft.Quantum.Chemistry.Fermion;
using Microsoft.Quantum.Chemistry.Magic;
using Newtonsoft.Json;
using Xunit;
using Microsoft.Quantum.Chemistry.JordanWigner;
using Microsoft.Quantum.Chemistry.OrbitalIntegrals;

namespace Microsoft.Quantum.Chemistry.Tests
{
    public class ChemistryEncodeMagicTests
    {
        public (ChemistryEncodeMagic, MockChannel) Init() =>
            (new ChemistryEncodeMagic(), new MockChannel());


        [Fact]
        public async Task EncodeBroombridge()
        {
            var (magic, channel) = Init();
            var fileName = "broombridge_v0.2.yaml";

            Assert.Equal("%chemistry.encode", magic.Name);

            using var reader = File.OpenText(fileName);
            var broombridge = BroombridgeSerializer.Deserialize(reader);
            var problemData = broombridge.First();
            var orbitalIntegralHamiltonian = problemData.OrbitalIntegralHamiltonian;
            var fermionHamiltonian = orbitalIntegralHamiltonian.ToFermionHamiltonian(IndexConvention.UpDown);

            var wavefunctionMagic = new WavefunctionMagic();
            var args = JsonConvert.SerializeObject(new WavefunctionMagic.Arguments
            {
                FileName = "broombridge_v0.2.yaml",
                WavefunctionLabel = "UCCSD |G>",
                IndexConvention = IndexConvention.HalfUp
            });
            var wavefunction = (FermionWavefunction<int>)((await wavefunctionMagic.Run(args, channel)).Output);

            args = JsonConvert.SerializeObject(new ChemistryEncodeMagic.Arguments
            {
                Hamiltonian = fermionHamiltonian,
                Wavefunction = wavefunction
            });

            var result = await magic.Run(args, channel);
            Assert.Equal(ExecuteStatus.Ok, result.Status);
            var data = result.Output as JordanWignerEncodingData;

            Assert.Equal(12, data.Item1);
        }
    }
}