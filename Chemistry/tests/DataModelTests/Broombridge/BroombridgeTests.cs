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

    public class BroombridgeVersionNumberTests
    {
        [Fact]
        public void DeserializeVersionNumbers()
        {
            Assert.Equal(Broombridge.Version.Number.v0_1, Broombridge.Deserialize.GetVersionNumber("Broombridge/broombridge_v0.1.yaml"));
            Assert.Equal(Broombridge.Version.Number.v0_2, Broombridge.Deserialize.GetVersionNumber("Broombridge/broombridge_v0.2.yaml"));
        }
    }

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
            Assert.Equal("0.2", broombridge.Format.Version);
            Assert.NotEqual("", broombridge.Format.Version);
        }

        [Fact]
        public void UnitaryCoupledCluster()
        {
            var state = broombridge.ProblemDescription.First().InitialStates.ElementAt(3);

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
        public void UpdateFrom_v0_1()
        {
            var filename = "Broombridge/broombridge_v0.1.yaml";
            var broombridge_v0_1 = Broombridge.Deserialize.v0_1(filename);
            var broombridge_v0_2 = Broombridge.Update.Data(broombridge_v0_1);

            Broombridge.Serialize.v0_2(broombridge_v0_2, "");

        }
    }

    public class BroombridgeTypedTests
    {
        static string filename = "Broombridge/broombridge_v0.2.yaml";
        static Broombridge.V0_2.Data broombridge = Broombridge.Deserialize.v0_2(filename);
        static Broombridge.V0_2.ProblemDescription broombridgeProblem = broombridge.ProblemDescription.First();
        
        [Fact]
        public void UnitaryCoupledCluster()
        {
            BroombridgeTyped broombridgeTyped = new BroombridgeTyped(broombridgeProblem);
            var state = broombridgeTyped.InitialStates["UCCSD |G>"];

            var targetTerm = new FermionTerm(
(long[])                (new[] { 1L, 1L, 0L, 0L }),
(SpinOrbital[])                (new[] { ((int)0, u:(Spin)Spin.u), ((int)1, u:(Spin)Spin.u), ((int)1, d:(Spin)Spin.d), ((int)3, d:(Spin)Spin.d) }).Select((Func<(int, Spin), SpinOrbital>)(((int, Spin) o) => (SpinOrbital)new SpinOrbital(o))).ToArray(),
(double)-0.5);

            var checkTerm = state.Superposition.ElementAt(2).term;

            Assert.Equal(targetTerm, checkTerm);
        }
    }
}