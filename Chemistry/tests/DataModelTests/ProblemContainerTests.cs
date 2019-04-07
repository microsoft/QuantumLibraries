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

   
    public class ProblemContainerTests
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