// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using System.Collections.Generic;
using System.IO;
using System.Linq;

using Microsoft.Quantum.Chemistry.Broombridge;

using Newtonsoft.Json;

using Xunit;

namespace Microsoft.Quantum.Chemistry.Tests
{
    public class BroombridgeVersionNumberTests
    {
        [Fact]
        public void DeserializeVersionNumbers()
        {
            Assert.Equal(VersionNumber.v0_1, Broombridge.Deserializers.GetVersionNumber("Broombridge/broombridge_v0.1.yaml"));
            Assert.Equal(VersionNumber.v0_2, Broombridge.Deserializers.GetVersionNumber("Broombridge/broombridge_v0.2.yaml"));
        }
    }

    public class Broombridgev0_1Tests
    {

        [Fact]
        public void LoadTest()
        {
            var filename = "Broombridge/broombridge_v0.1.yaml";


            var broombridge = Broombridge.Deserializers.Deserialize<V0_1.Data>(filename);

            Assert.Equal("0.1", broombridge.Format.Version);
        }
    }

    public class Broombridgev0_2Tests
    {
        static string filename = "Broombridge/broombridge_v0.2.yaml";
        static V0_2.Data broombridge = Deserializers.Deserialize<V0_2.Data>(filename);


        [Fact]
        public void Version()
        {
            Assert.Equal("0.2", broombridge.Format.Version);
            Assert.NotEqual("", broombridge.Format.Version);
        }

        [Fact]
        public void UnitaryCoupledCluster()
        {
            var state = broombridge.ProblemDescriptions.Single().InitialStates.ElementAt(3);

            Assert.Equal("UCCSD |G>", state.Label);

            Assert.Equal("unitary_coupled_cluster", state.Method);

            //Assert.Equal("|G>", state.ClusterOperator.Reference);
            Assert.Equal(new[] {"1.0", "(1a)+", "(2a)+", "(3a)+", "(4a)+", "(5a)+", "(6a)+", "(1b)+", "(2b)+", "(3b)+", "(4b)+", "(5b)+", "(6b)+", "|vacuum>"}, state.ClusterOperator.Reference);

            var oneBodyAmplitude = state.ClusterOperator.OneBodyAmplitudes;
            List<List<string>> oneBodyAmplitudeTruth =
                new[] {
                    new[] {"0.1", "(1a)+", "(2a)"}.ToList(),
                    new[] {"-0.2", "(1b)+", "(2a)"}.ToList()
                }.ToList();
            Assert.Equal(oneBodyAmplitudeTruth, oneBodyAmplitude);

            var twoBodyAmplitude = state.ClusterOperator.TwoBodyAmplitudes;
            List<List<string>> twoBodyAmplitudeTruth =
                new[] {
                    new[] {"-0.5", "(1a)+", "(2a)+", "(2b)", "(4b)"}.ToList(),
                    new[] {"0.5", "(1a)+", "(3a)+", "(2a)", "(4b)"}.ToList()
                }.ToList();
            Assert.Equal(twoBodyAmplitudeTruth, twoBodyAmplitude);

        }

        [Fact]
        public void CheckNullUnitaryCoupledCluster()
        {
            var broombridge_internal = Deserializers.DeserializeBroombridge(filename).ProblemDescriptions.Single();

            Assert.Contains("UCCSD nullTwo", broombridge_internal.Wavefunctions.Keys);
            Assert.Contains("UCCSD nullOne", broombridge_internal.Wavefunctions.Keys);
        }

        [Fact]
        public void UpdateFrom_v0_1()
        {
            var filename = "Broombridge/broombridge_v0.1.yaml";
            var broombridge_v0_1 = Deserializers.Deserialize<V0_1.Data>(filename);
            var broombridge_v0_2 = DataStructures.Update(broombridge_v0_1);

            Broombridge.Serializers.SerializeBroombridgev0_2(broombridge_v0_2, "");

        }

        [Fact]
        public void JsonEncoding()
        {
            var filename = "Broombridge/broombridge_v0.2.yaml";
            var original = Deserializers.DeserializeBroombridge(filename).Raw;
            
            var json = JsonConvert.SerializeObject(original);
            File.WriteAllText("original.json", json);

            var serialized = JsonConvert.DeserializeObject<V0_2.Data>(json);
                File.WriteAllText("serialized.json", JsonConvert.SerializeObject(serialized));

            Assert.Equal(original.Format, serialized.Format);
            Assert.Equal(original.Bibliography.Count, serialized.Bibliography.Count);
            Assert.Equal(original.ProblemDescriptions.Count, serialized.ProblemDescriptions.Count);
            Assert.Equal(original.Generator.Source, serialized.Generator.Source);
            Assert.Equal(original.Schema, serialized.Schema);
        }
    }
}