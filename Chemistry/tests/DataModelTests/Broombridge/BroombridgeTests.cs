// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using Xunit;
using Microsoft.Quantum.Chemistry;
using Microsoft.Quantum.Simulation.Core;
using System.Text.RegularExpressions;
using System;
using System.Linq;
using System.Collections.Generic;
using System.IO;

using YamlDotNet.Core;
using YamlDotNet.Serialization;

namespace Microsoft.Quantum.Chemistry.Tests
{
    using static FermionTermType.Common;
    using FermionTerm = FermionTerm;
    using FermionTermType = FermionTermType;
    using SpinOrbital = SpinOrbital;
    
    public class Broombridgev0_1Tests
    {

        [Fact]
        public void LoadTest()
        {
            var filename = "Broombridge/broombridge_v0.1.yaml";
            

            var broombridge = Broombridge.Deserialize.v0_1(filename);

            Assert.Equal("0.1", broombridge.Format.Version);
        }
    }

    public class Broombridgev0_2Tests
    {
        static string filename = "Broombridge/broombridge_v0.2.yaml";
        static Broombridge.V0_2.Data broombridge = Broombridge.Deserialize.v0_2(filename);

        [Fact]
        public void Version()
        {
            Assert.Equal("0.2", broombridge.Version);
            Assert.NotEqual("", broombridge.Version);
        }

        [Fact]
        public void UnitaryCoupledCluster()
        {
            var state = broombridge.ProblemDescription.First().InitialStates.ElementAt(3);
            
            Assert.Equal("UCCSD |G>", state.Label);

            Assert.Equal("unitary_coupled_cluster", state.Method);

            Assert.Equal("|G>", state.ClusterOperator.Reference);

            var oneBodyAmplitude = state.ClusterOperator.OneBodyAmplitudes;
            List<List<string>> oneBodyAmplitudeTruth =
                new[] {
                    new[] { "1a" , "2a", "0.1"}.ToList(),
                    new[] { "1b" , "2a", "-0.2"}.ToList()
                }.ToList();
            Assert.Equal(oneBodyAmplitudeTruth, oneBodyAmplitude);

            var twoBodyAmplitude = state.ClusterOperator.TwoBodyAmplitudes;
            List<List<string>> twoBodyAmplitudeTruth =
                new[] {
                    new[] { "1a" , "2a", "2b", "4b", "-0.5"}.ToList(),
                    new[] { "1a" , "3a", "2a", "4b", "0.5"}.ToList()
                }.ToList();
            Assert.Equal(twoBodyAmplitudeTruth, twoBodyAmplitude);

        }
    }


}