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

namespace Microsoft.Quantum.Chemistry.Tests
{
    public class FermionHamiltonianLoadMagicTests
    {
        public (FermionHamiltonianLoadMagic, MockChannel) Init() =>
            (new FermionHamiltonianLoadMagic(), new MockChannel());

        [Fact]
        public void LoadNoFile()
        {
            var (magic, channel) = Init();

            Assert.Equal("%fh_load", magic.Name);
            var result = magic.Run("", channel);

            Assert.Equal(ExecuteStatus.Error, result.Status);
        }


        [Fact]
        public void LoadInvalidFile()
        {
            var (magic, channel) = Init();
            var filename = "Broombridge/foo_bar.yaml";

            Assert.Throws<FileNotFoundException>(() => magic.Run(filename, channel));
        }


        [Fact]
        public void LoadBroombridgeFile()
        {
            var (magic, channel) = Init();
            var filename = "Broombridge/broombridge_v0.2.yaml";

            var result = magic.Run(filename, channel);
            var hamiltonian = result.Output as FermionHamiltonian;
            Assert.Equal(ExecuteStatus.Ok, result.Status);
            Assert.Equal(12, hamiltonian.SystemIndices.Count);
            Assert.Equal(6, hamiltonian.Terms.Count);
            Assert.Equal(64.5730917943, hamiltonian.Norm());
        }
    }

    public class FermionHamiltonianCreateMagicTests
    {
        public (FermionHamiltonianCreateMagic, MockChannel) Init() =>
            (new FermionHamiltonianCreateMagic(), new MockChannel());

        [Fact]
        public void CreateHamiltonian()
        {
            var (magic, channel) = Init();

            Assert.Equal("%fh_create", magic.Name);
            var result = magic.Run("", channel);
            var hamiltonian = result.Output as FermionHamiltonian;
            Assert.Equal(ExecuteStatus.Ok, result.Status);
            Assert.Empty(hamiltonian.SystemIndices);
            Assert.Empty(hamiltonian.Terms);
        }
    }

    public class FermionHamiltonianAddTermsMagicTests
    {
        public (FermionHamiltonianCreateMagic, FermionHamiltonianAddTermsMagic, MockChannel) Init() =>
            (new FermionHamiltonianCreateMagic(), new FermionHamiltonianAddTermsMagic(), new MockChannel());

        [Fact]
        public void AddHamiltonianTerms()
        {
            var (createMagic, magic, channel) = Init();
            Assert.Equal("%fh_add_terms", magic.Name);

            var original = createMagic.Run("", channel).Output as FermionHamiltonian;
            var fermionTerms = new List<(HermitianFermionTerm, double)>
                {
                    (new HermitianFermionTerm(new int[] {}.ToLadderSequence()), 10.0),
                    (new HermitianFermionTerm(new[] {0,0}.ToLadderSequence()), 1.0),
                    (new HermitianFermionTerm(new[] {1,1}.ToLadderSequence()), 1.0),
                    (new HermitianFermionTerm(new[] {2,2}.ToLadderSequence()), 1.0),
                    (new HermitianFermionTerm(new[] {0,2}.ToLadderSequence()), 1.0),
                    (new HermitianFermionTerm(new[] {1,3}.ToLadderSequence()), 1.0),
                    (new HermitianFermionTerm(new[] {2,6}.ToLadderSequence()), 1.0),
                    (new HermitianFermionTerm(new[] {0,2,2,0}.ToLadderSequence()), 1.0),
                    (new HermitianFermionTerm(new[] {1,3,3,1}.ToLadderSequence()), 1.0),
                    (new HermitianFermionTerm(new[] {2,6,6,2}.ToLadderSequence()), 1.0),
                    (new HermitianFermionTerm(new[] {0,2,2,1}.ToLadderSequence()), 1.0),
                    (new HermitianFermionTerm(new[] {1,3,3,2}.ToLadderSequence()), 1.0),
                    (new HermitianFermionTerm(new[] {2,6,6,5}.ToLadderSequence()), 1.0),
                    (new HermitianFermionTerm(new[] {0,2,4,3}.ToLadderSequence()), 1.0),
                    (new HermitianFermionTerm(new[] {1,4,3,2}.ToLadderSequence()), 1.0),
                    (new HermitianFermionTerm(new[] {2,4,5,3}.ToLadderSequence()), 1.0)
                };

            var args = JsonConvert.SerializeObject(new FermionHamiltonianAddTermsMagic.Arguments
            {
                hamiltonian = original,
                fermionTerms = fermionTerms
            });

            var result = magic.Run(args, channel);
            var hamiltonian = result.Output as FermionHamiltonian;
            Assert.Equal(ExecuteStatus.Ok, result.Status);
            Assert.Equal(16, hamiltonian.CountTerms());
        }
    }
}