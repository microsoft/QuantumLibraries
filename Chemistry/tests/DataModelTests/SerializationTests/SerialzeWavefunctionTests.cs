// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Diagnostics;
using Microsoft.Quantum.Chemistry.Broombridge;
using Microsoft.Quantum.Chemistry.Fermion;
using Newtonsoft.Json;
using Microsoft.Quantum.Chemistry.LadderOperators;
using Xunit;
using Microsoft.Quantum.Chemistry.LadderOperators;
using Microsoft.Quantum.Chemistry.JordanWigner;
using Microsoft.Quantum.Chemistry.OrbitalIntegrals;

namespace Microsoft.Quantum.Chemistry.Tests
{
    public class SerializeWavefunctionTests
    {
        public FermionWavefunction<int> MakeWavefunction(string wavefunctionLabel = "UCCSD |G>")
        {
            var data = Broombridge.Deserializers.DeserializeBroombridge(
                Path.Join("Broombridge", "broombridge_v0.2.yaml"
            ));
            return data.ProblemDescriptions.First().Wavefunctions[wavefunctionLabel].ToIndexing(IndexConvention.UpDown);
        }

        [Fact]
        public void LoadWavefunctionTest()
        {
            var wavefunction = MakeWavefunction("UCCSD |G>");
        }

        [Fact]
        public void SerializeUCCWavefunctionTest()
        {
            var wavefunction = MakeWavefunction("UCCSD |G>");

            string json = JsonConvert.SerializeObject(wavefunction, Formatting.Indented);

            Debug.WriteLine(json);

            var deserializedWavefunction = JsonConvert.DeserializeObject<FermionWavefunction<int>>(json);
        }

        [Fact]
        public void SerializeMCFWavefunctionTest()
        {
            var wavefunction = MakeWavefunction("|E>");

            string json = JsonConvert.SerializeObject(wavefunction, Formatting.Indented);

            Debug.WriteLine(json);

            var deserializedWavefunction = JsonConvert.DeserializeObject<FermionWavefunction<int>>(json);
        }


    }
}