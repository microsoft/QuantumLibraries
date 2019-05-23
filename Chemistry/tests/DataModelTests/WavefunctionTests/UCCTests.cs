// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using Xunit;
using Microsoft.Quantum.Chemistry;
using Microsoft.Quantum.Simulation.Core;

using System.Text.RegularExpressions;
using System;
using System.Linq;
using System.Collections.Generic;


using Microsoft.Quantum.Chemistry.Fermion;
using Microsoft.Quantum.Chemistry.OrbitalIntegrals;
using Microsoft.Quantum.Chemistry.LadderOperators;

namespace Microsoft.Quantum.Chemistry.Tests
{
    using HermitianFermionTerm = HermitianFermionTerm;
    using static TermType.Fermion;
    using static RaisingLowering;

    public class WavefunctionTests
    {
        [Fact]
        public void MakeUCCSDExcitationsTest0()
        {
            var reference = new SingleCFWavefunction<SpinOrbital>(new[] { (1, Spin.u) }.Select(o => new SpinOrbital(o)));
            var state = reference.AddAllUCCSDSingletExcitations(3);

            var excitations = new[]
            {
                new []{ (0, Spin.u), (1,Spin.u)},
                new []{ (2, Spin.u), (1,Spin.u)}
            }. Select(o => new IndexOrderedSequence<SpinOrbital>(o.ToLadderSequence()));


            Assert.True(state.Excitations.Keys.All(excitations.Contains));
            Assert.True(state.Excitations.Count() == excitations.Count());
        }

        [Fact]
        public void MakeUCCSDExcitationsTest1()
        {
            var reference = new SingleCFWavefunction<SpinOrbital>(new[] { (1, Spin.d), (2, Spin.d) }.Select(o => new SpinOrbital(o)));
            var state = reference.AddAllUCCSDSingletExcitations(3);

            var excitations = new[]
            {
                new []{ (0, Spin.d), (1,Spin.d)},
                new []{ (0, Spin.d), (2,Spin.d)}
            }.Select(o => new IndexOrderedSequence<SpinOrbital>(o.ToLadderSequence()));


            Assert.True(state.Excitations.Keys.All(excitations.Contains));
            Assert.True(state.Excitations.Count() == excitations.Count());
        }

        [Fact]
        public void MakeUCCSDExcitationsTest2()
        {
            var reference = new SingleCFWavefunction<SpinOrbital>(new[] { (1, Spin.u), (2, Spin.d) }.Select(o => new SpinOrbital(o)));
            var state = reference.AddAllUCCSDSingletExcitations(3);

            var excitations = new[]
            {
                new []{ (0, Spin.u), (1,Spin.u)},
                new []{ (2, Spin.u), (1,Spin.u)},
                new []{ (0, Spin.d), (2,Spin.d)},
                new []{ (1, Spin.d), (2,Spin.d)},
                new []{ (0, Spin.u), (0, Spin.d), (2, Spin.d), (1,Spin.u)},
                new []{ (0, Spin.u), (1, Spin.d), (2, Spin.d), (1,Spin.u)},
                new []{ (0, Spin.d), (2, Spin.u), (2, Spin.d), (1,Spin.u)},
                new []{ (1, Spin.d), (2, Spin.u), (2, Spin.d), (1,Spin.u)}
            }.Select(o => new IndexOrderedSequence<SpinOrbital>(o.ToLadderSequence()));

            Assert.True(state.Excitations.Keys.All(excitations.Contains));
            Assert.True(state.Excitations.Count() == excitations.Count());
        }
    }

}