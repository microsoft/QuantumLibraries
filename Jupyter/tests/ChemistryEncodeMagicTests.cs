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
using Microsoft.Quantum.Chemistry.LadderOperators;
using Microsoft.Quantum.Chemistry.JordanWigner;
using Microsoft.Quantum.Chemistry.OrbitalIntegrals;

namespace Microsoft.Quantum.Chemistry.Tests
{
    public class ChemistryEncodeMagicTests
    {
        public (ChemistryEncodeMagic, MockChannel) Init() =>
            (new ChemistryEncodeMagic(), new MockChannel());


        [Fact]
        public void EncodeBroombridge()
        {
            var (magic, channel) = Init();
            var fileName = "broombridge_v0.2.yaml";

            Assert.Equal("%chemistry.encode", magic.Name);

            var broombridge = Deserializers.DeserializeBroombridge(fileName);
            var problemData = broombridge.ProblemDescriptions.First();
            var orbitalIntegralHamiltonian = problemData.OrbitalIntegralHamiltonian;
            var fermionHamiltonian = orbitalIntegralHamiltonian.ToFermionHamiltonian(IndexConvention.UpDown);

            var wavefunctionMagic = new WavefunctionMagic();
            var args = JsonConvert.SerializeObject(new WavefunctionMagic.Arguments
            {
                FileName = "broombridge_v0.2.yaml",
                WavefunctionLabel = "UCCSD |G>",
                IndexConvention = IndexConvention.HalfUp
            });
            var wavefunction = (FermionWavefunction<int>)wavefunctionMagic.Run(args, channel).Output;

            args = JsonConvert.SerializeObject(new ChemistryEncodeMagic.Arguments
            {
                hamiltonian = fermionHamiltonian,
                wavefunction = wavefunction
            });

            var result = magic.Run(args, channel);
            Assert.Equal(ExecuteStatus.Ok, result.Status);
            var data = result.Output as JordanWignerEncodingData;

            Assert.Equal(12, data.Item1);
        }
    }
}